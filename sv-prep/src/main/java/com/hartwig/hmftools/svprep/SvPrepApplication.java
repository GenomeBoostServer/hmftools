package com.hartwig.hmftools.svprep;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.TaskExecutor.runThreadTasks;
import static com.hartwig.hmftools.svprep.SvCommon.APP_NAME;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.svprep.reads.PartitionTask;
import com.hartwig.hmftools.svprep.reads.PartitionThread;

import org.jetbrains.annotations.NotNull;

public class SvPrepApplication
{
    private final SvConfig mConfig;
    private final ResultsWriter mWriter;
    private final SpanningReadCache mSpanningReadCache;
    private final ExistingJunctionCache mExistingJunctionCache;

    public SvPrepApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new SvConfig(configBuilder);
        mWriter = new ResultsWriter(mConfig);
        mSpanningReadCache = new SpanningReadCache(mConfig);
        mExistingJunctionCache = new ExistingJunctionCache();
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        SV_LOGGER.info("running SvPrep for sample({})", mConfig.SampleId);

        long startTimeMs = System.currentTimeMillis();

        if(mConfig.CalcFragmentLength)
            calcFragmentDistribution();

        mExistingJunctionCache.loadJunctions(mConfig.ExistingJunctionFile);

        CombinedStats combinedStats = new CombinedStats();

        if(mConfig.GlobalQueue)
        {
            processGlobalQueue(combinedStats);
        }
        else
        {
            processPerChromosome(combinedStats);
        }

        if(mConfig.UseCacheBam)
        {
            mSpanningReadCache.reset(); // clear data before candidate assignment
            mSpanningReadCache.candidateBamWriter().assignCandidateReads(mWriter);
        }

        mWriter.close();

        long timeTakenMs = System.currentTimeMillis() - startTimeMs;

        if(mConfig.PerfDebug && (combinedStats.ReadStats.TotalReads > 10000 || timeTakenMs > 10000))
        {
            SV_LOGGER.info("final stats: {}", combinedStats.ReadStats.toString());

            if(SV_LOGGER.isDebugEnabled())
                combinedStats.PerfCounters.forEach(x -> x.logIntervalStats());
            else
                combinedStats.PerfCounters.forEach(x -> x.logStats());
        }

        SV_LOGGER.info("SvPrep complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    // Legacy path — use -per_chromosome flag to activate (creates a thread pool per chromosome).
    private void processPerChromosome(final CombinedStats combinedStats)
    {
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                continue;

            SV_LOGGER.info("processing chromosome({})", chromosomeStr);

            ChromosomeTask chromosomeTask = new ChromosomeTask(chromosomeStr, mConfig, mSpanningReadCache, mExistingJunctionCache, mWriter);
            chromosomeTask.process();
            combinedStats.addPartitionStats(chromosomeTask.combinedStats().ReadStats);

            if(combinedStats.PerfCounters.isEmpty())
            {
                combinedStats.PerfCounters.addAll(chromosomeTask.combinedStats().PerfCounters);
            }
            else
            {
                for(int i = 0; i < chromosomeTask.combinedStats().PerfCounters.size(); ++i)
                {
                    combinedStats.PerfCounters.get(i).merge(chromosomeTask.combinedStats().PerfCounters.get((i)));
                }
            }
        }
    }

    // Drain all partitions from a single global queue using a single thread pool.
    private void processGlobalQueue(final CombinedStats combinedStats)
    {
        List<ChrBaseRegion> allRegions = new ArrayList<>();
        int chromosomeCount = 0;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                continue;

            allRegions.addAll(ChromosomeTask.buildPartitions(chromosomeStr, mConfig));
            ++chromosomeCount;
        }

        // Queue partitions in reference order.
        Queue<PartitionTask> partitions = new ConcurrentLinkedQueue<>();

        int taskId = 0;
        for(ChrBaseRegion region : allRegions)
        {
            partitions.add(new PartitionTask(region, taskId++));
        }

        int workerCount = min(partitions.size(), mConfig.Threads);

        SV_LOGGER.info("global partition queue: partitions({}) chromosomes({}) workers({})",
                partitions.size(), chromosomeCount, workerCount);

        List<Thread> workers = new ArrayList<>();

        for(int i = 0; i < workerCount; ++i)
        {
            workers.add(new PartitionThread(
                    "ALL", mConfig, partitions, mSpanningReadCache,
                    mExistingJunctionCache, mWriter, combinedStats));
        }

        if(!runThreadTasks(workers))
            System.exit(1);

        SV_LOGGER.info("global queue complete, stats: {}", combinedStats.ReadStats.toString());

        mSpanningReadCache.logStats();
    }

    private void calcFragmentDistribution()
    {
        FragmentSizeDistribution fragSizeDistribution = new FragmentSizeDistribution(mConfig);
        fragSizeDistribution.run();

        final int[] fragmentDistLengths = fragSizeDistribution.calculatePercentileLengths();

        SV_LOGGER.info("fragment length distribution percentile values(min={} max={})",
                fragmentDistLengths[0], fragmentDistLengths[1]);

        if(fragmentDistLengths[0] > 0 && fragmentDistLengths[1] > 0)
            mConfig.ReadFiltering.config().setFragmentLengths(fragmentDistLengths[0], fragmentDistLengths[1]);

        if(fragSizeDistribution.maxReadLength() > 0)
            mConfig.ReadLength = fragSizeDistribution.maxReadLength();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SvConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SvPrepApplication svPrep = new SvPrepApplication(configBuilder);
        svPrep.run();
    }
}
