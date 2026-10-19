package com.hartwig.hmftools.svprep.depth;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.svprep.SvCommon.APP_NAME;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.ExcludedRegions;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.svprep.reads.PartitionTask;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class DepthAnnotator
{
    private final DepthConfig mConfig;
    private final Map<String,Integer> mSampleVcfGenotypeIds;
    private final List<ChrBaseRegion> mExcludedRegions;

    public DepthAnnotator(final ConfigBuilder configBuilder)
    {
        mConfig = new DepthConfig(configBuilder);
        mSampleVcfGenotypeIds = Maps.newHashMap();
        mExcludedRegions = ExcludedRegions.getPolyGRegions(mConfig.RefGenVersion);
    }

    public void run()
    {
        if(mConfig.InputVcf == null || !Files.exists(Paths.get(mConfig.InputVcf)) || mConfig.OutputVcf == null)
        {
            SV_LOGGER.error("missing VCF config or file");
            System.exit(1);
        }

        if(mConfig.BamFiles.size() != mConfig.Samples.size())
        {
            SV_LOGGER.error("inconsistent samples and BAM files");
            System.exit(1);
        }

        SV_LOGGER.info("SvPrep depth annotation for samples: {}", mConfig.Samples);

        long startTimeMs = System.currentTimeMillis();

        VcfFileReader reader = new VcfFileReader(mConfig.InputVcf);
        VCFHeader vcfHeader = reader.vcfHeader();

        if(!establishGenotypeIds(vcfHeader))
            System.exit(1);

        reader.close(); // each worker opens its own reader for concurrent region queries

        ensureVcfIndex(vcfHeader); // region queries below require a tabix index

        // partition the genome into windows, in coordinate order, so only a bounded number of
        // windows is ever resident in memory at once instead of the entire VCF
        List<ChrBaseRegion> windows = buildWindows();

        Queue<PartitionTask> partitions = new ConcurrentLinkedQueue<>();
        int taskId = 0;

        for(ChrBaseRegion region : windows)
        {
            partitions.add(new PartitionTask(region, taskId++));
        }

        int regionCount = partitions.size();

        int threadCount = max(mConfig.Threads, 1);

        // cap of windows held in memory at once (>= threadCount, plus a little read-ahead);
        // kept tight because each window holds its variants and (while slicing) its BAM reads
        int maxPendingRegions = max(threadCount + 4, 8);

        DepthWriter depthWriter = new DepthWriter(mConfig, vcfHeader, maxPendingRegions);

        if(regionCount == 0)
        {
            // still write a header-only VCF (as the full-load version did) so downstream has input
            SV_LOGGER.warn("no regions to process");
            depthWriter.close();
            return;
        }

        SV_LOGGER.info("depth annotation: {} regions, threads({}) partitionSize({}) maxPending({})",
                regionCount, threadCount, mConfig.PartitionSize, maxPendingRegions);

        List<DepthPartitionThread> workers = Lists.newArrayList();

        for(int i = 0; i < min(regionCount, threadCount); ++i)
        {
            workers.add(new DepthPartitionThread(
                    mConfig, partitions, mSampleVcfGenotypeIds, mExcludedRegions, depthWriter));
        }

        if(!runThreadTasks(new ArrayList<>(workers)))
            System.exit(1);

        depthWriter.close();

        SV_LOGGER.info("SvPrep depth annotation complete, mins({})", runTimeMinsStr(startTimeMs));

        if(!workers.isEmpty())
        {
            PerformanceCounter perfCounter = workers.get(0).getPerfCounter();

            for(int i = 1; i < workers.size(); ++i)
            {
                perfCounter.merge(workers.get(i).getPerfCounter());
            }

            perfCounter.logStats();
        }
    }

    private List<ChrBaseRegion> buildWindows()
    {
        // with specific regions set (a testing/debug path) fall back to fixed-size partitioning
        if(!mConfig.SpecificRegions.isEmpty())
        {
            List<ChrBaseRegion> regions = Lists.newArrayList();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());
                regions.addAll(partitionChromosome(chrStr, mConfig.RefGenVersion,
                        mConfig.SpecificRegions, mConfig.PartitionSize));
            }

            return regions;
        }

        return buildGapAlignedWindows();
    }

    // Cut window boundaries only at gaps larger than ProximityDistance, so no proximity group (the
    // unit DepthTask slices) is ever split and output matches the full-load run exactly. A window
    // closes once it hits PartitionSize and the next gap appears. Positions from a light bgzf scan.
    private List<ChrBaseRegion> buildGapAlignedWindows()
    {
        // rank of each chromosome in HumanChromosome order (the order the full-load version wrote)
        Map<String,Integer> chromosomeRank = Maps.newHashMap();
        for(HumanChromosome chromosome : HumanChromosome.values())
            chromosomeRank.put(mConfig.RefGenVersion.versionedChromosome(chromosome.toString()),
                    chromosome.ordinal());

        Set<String> seenChromosomes = new HashSet<>();
        List<ChrBaseRegion> windows = Lists.newArrayList();

        String currentChromosome = null;
        int windowStart = -1;
        int prevPosition = -1;

        try(BufferedReader reader = new BufferedReader(
                new InputStreamReader(new BlockCompressedInputStream(new File(mConfig.InputVcf)))))
        {
            String line;
            while((line = reader.readLine()) != null)
            {
                if(line.isEmpty() || line.charAt(0) == '#')
                    continue;

                int tab1 = line.indexOf('\t');
                int tab2 = line.indexOf('\t', tab1 + 1);
                String chromosome = line.substring(0, tab1);

                if(!chromosomeRank.containsKey(chromosome))
                    continue;

                int position = Integer.parseInt(line.substring(tab1 + 1, tab2));

                if(!chromosome.equals(currentChromosome))
                {
                    if(currentChromosome != null)
                        windows.add(
                                new ChrBaseRegion(currentChromosome, windowStart, prevPosition));

                    // a contig split into separate blocks would produce overlapping windows and
                    // duplicate records
                    if(!seenChromosomes.add(chromosome))
                    {
                        SV_LOGGER.error("input VCF({}) not sorted: {} in separate blocks",
                                mConfig.InputVcf, chromosome);
                        System.exit(1);
                    }

                    currentChromosome = chromosome;
                    windowStart = position;
                    prevPosition = position;
                    continue;
                }

                if(prevPosition - windowStart + 1 >= mConfig.PartitionSize
                && position - prevPosition > mConfig.ProximityDistance)
                {
                    windows.add(new ChrBaseRegion(currentChromosome, windowStart, prevPosition));
                    windowStart = position;
                }

                prevPosition = position;
            }

            if(currentChromosome != null)
                windows.add(new ChrBaseRegion(currentChromosome, windowStart, prevPosition));
        }
        catch(Exception e)
        {
            SV_LOGGER.error("failed to pre-scan VCF({}) for window boundaries: {}",
                    mConfig.InputVcf, e.toString());
            System.exit(1);
        }

        // order windows by chromosome (HumanChromosome order) regardless of input contig order;
        // stable sort keeps each chromosome's windows in coordinate order (no-op if sorted)
        windows.sort(Comparator.comparingInt(window -> chromosomeRank.get(window.Chromosome)));

        return windows;
    }

    // region (tabix) queries need an index next to the bgzipped input; create it if missing, and
    // rebuild it if older than the VCF (a stale index would silently miss variants)
    private void ensureVcfIndex(final VCFHeader header)
    {
        File vcfFile = new File(mConfig.InputVcf);
        File indexFile = new File(mConfig.InputVcf + TABIX_EXTENSION);

        if(indexFile.exists())
        {
            if(indexFile.lastModified() >= vcfFile.lastModified())
                return;

            SV_LOGGER.warn("tabix index({}) is older than input VCF, rebuilding", indexFile);
        }

        if(!mConfig.InputVcf.endsWith(".gz"))
        {
            SV_LOGGER.error("input VCF({}) must be bgzipped and tabix-indexed for region queries",
                    mConfig.InputVcf);
            System.exit(1);
        }

        SV_LOGGER.info("creating tabix index for input VCF({})", mConfig.InputVcf);

        try
        {
            Path vcfPath = Paths.get(mConfig.InputVcf);

            TabixIndex index = IndexFactory.createTabixIndex(
                    vcfPath.toFile(), new VCFCodec(), TabixFormat.VCF,
                    header.getSequenceDictionary());

            index.writeBasedOnFeaturePath(vcfPath);
        }
        catch(Exception e)
        {
            SV_LOGGER.error("failed to create tabix index for input VCF({}): {}",
                    mConfig.InputVcf, e.toString());
            System.exit(1);
        }
    }

    private static final String TABIX_EXTENSION = ".tbi";

    private boolean establishGenotypeIds(final VCFHeader header)
    {
        List<String> vcfSampleNames = header.getGenotypeSamples();

        for(int s = 0; s < mConfig.Samples.size(); ++s)
        {
            String sampleId = mConfig.Samples.get(s);
            boolean found = false;

            for(int i = 0; i < vcfSampleNames.size(); ++i)
            {
                String vcfSampleName = vcfSampleNames.get(i);

                if(sampleId.equals(vcfSampleName))
                {
                    mSampleVcfGenotypeIds.put(sampleId, i);
                    found = true;
                    break;
                }
            }

            if(!found)
            {
                SV_LOGGER.error("sample({}) not found in VCF genotype names", sampleId);
                return false;
            }
        }

        return true;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        DepthConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        DepthAnnotator bamSvSlicer = new DepthAnnotator(configBuilder);
        bamSvSlicer.run();
    }
}
