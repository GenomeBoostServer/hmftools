package com.hartwig.hmftools.svprep.depth;

import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Queue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.svprep.reads.PartitionTask;

import htsjdk.samtools.SamReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

// one worker per thread: pulls windows from the shared queue, loads only that window's variants via
// a tabix query, annotates them (DepthTask) and writes through the ordered writer. Each worker owns
// its own readers because htsjdk query iterators are not safe to share across threads.
public class DepthPartitionThread extends Thread
{
    private final DepthConfig mConfig;
    private final Queue<PartitionTask> mPartitions;
    private final Map<String,Integer> mSampleVcfGenotypeIds;
    private final List<ChrBaseRegion> mExcludedRegions;
    private final DepthWriter mDepthWriter;

    private final VcfFileReader mVcfReader;
    private final List<SamReader> mSamReaders;

    private final PerformanceCounter mPerfCounter;

    public DepthPartitionThread(
            final DepthConfig config, final Queue<PartitionTask> partitions,
            final Map<String,Integer> sampleVcfGenotypeIds,
            final List<ChrBaseRegion> excludedRegions, final DepthWriter depthWriter)
    {
        mConfig = config;
        mPartitions = partitions;
        mSampleVcfGenotypeIds = sampleVcfGenotypeIds;
        mExcludedRegions = excludedRegions;
        mDepthWriter = depthWriter;

        mVcfReader = new VcfFileReader(mConfig.InputVcf, true);
        mSamReaders = DepthTask.openSamReaders(mConfig);

        mPerfCounter = new PerformanceCounter("Slice");

        start();
    }

    public PerformanceCounter getPerfCounter() { return mPerfCounter; }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                // acquire the permit before taking a task, so the next region to write is held
                // by a worker and the writer can flush (taking a task first could deadlock)
                mDepthWriter.acquireRegionPermit();

                PartitionTask partition = mPartitions.poll();

                if(partition == null)
                {
                    mDepthWriter.releaseRegionPermit();
                    SV_LOGGER.trace("all depth regions complete");
                    break;
                }

                processRegion(partition);
            }
            catch(Throwable e)
            {
                SV_LOGGER.error("depth thread execution error: {}", e.toString());
                e.printStackTrace();
                System.exit(1);
            }
        }

        mVcfReader.close();

        for(SamReader samReader : mSamReaders)
        {
            try
            {
                samReader.close();
            }
            catch(IOException e)
            {
                SV_LOGGER.error("failed to close bam file: {}", e.toString());
            }
        }
    }

    private void processRegion(final PartitionTask partition)
    {
        ChrBaseRegion region = partition.Region;

        List<VariantContext> variants = Lists.newArrayList();

        CloseableTribbleIterator<VariantContext> iterator = mVcfReader.regionIterator(region);

        // regionIterator logs and returns null on a failed query; treating that as an empty window
        // would silently drop its variants from the output, so fail the run instead
        if(iterator == null)
            throw new IllegalStateException("failed to query input VCF for region " + region);

        for(VariantContext variantContext : iterator)
        {
            if(excludeVariant(variantContext))
                continue;

            // a tabix query returns any variant overlapping the window, so assign it to the single
            // window whose range contains its start - prevents double-counting at window boundaries
            if(variantContext.getStart() < region.start() ||
                    variantContext.getStart() > region.end())
                continue;

            // create a new instance to set depth values against (matches the full-load behaviour)
            VariantContext newVariant = new VariantContextBuilder(variantContext)
                    .genotypes(variantContext.getGenotypes())
                    .filters(variantContext.getFilters())
                    .make();

            variants.add(newVariant);
        }

        iterator.close();

        if(variants.isEmpty())
        {
            mDepthWriter.writeRegion(partition.TaskId, variants);
            return;
        }

        DepthTask depthTask = new DepthTask(
                region.Chromosome, mConfig, mSampleVcfGenotypeIds, mSamReaders);
        depthTask.addVariants(variants);
        depthTask.call();

        mPerfCounter.merge(depthTask.getPerfCounter());

        mDepthWriter.writeRegion(partition.TaskId, depthTask.variants());
    }

    private boolean excludeVariant(final VariantContext variant)
    {
        if(mExcludedRegions.stream().anyMatch(x -> x.containsPosition(variant.getStart())))
            return true;

        if(mConfig.SpecificRegions.isEmpty())
            return false;

        return mConfig.SpecificRegions.stream()
                .noneMatch(x -> x.containsPosition(variant.getContig(), variant.getStart()));
    }
}
