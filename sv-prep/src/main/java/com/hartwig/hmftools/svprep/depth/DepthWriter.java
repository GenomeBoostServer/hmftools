package com.hartwig.hmftools.svprep.depth;

import static com.hartwig.hmftools.common.sv.SvVcfTags.ALLELE_FRACTION;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ALLELE_FRACTION_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR_DESC;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.util.List;
import java.util.TreeMap;
import java.util.concurrent.Semaphore;

import com.google.common.collect.Maps;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

// writes variants to the output VCF in coordinate order even though regions finish out of order:
// completed regions wait in an index-keyed reorder buffer until their turn. A semaphore bounds the
// regions resident in memory so a single slow region cannot re-grow the heap.
public class DepthWriter
{
    private final VariantContextWriter mWriter;

    private final Semaphore mPendingRegions;
    private final TreeMap<Integer,List<VariantContext>> mPendingVariants;
    private int mNextWriteIndex;
    private final Object mWriteLock;

    public DepthWriter(final DepthConfig config, final VCFHeader header, int maxPendingRegions)
    {
        mWriter = new VariantContextWriterBuilder()
                .setReferenceDictionary(header.getSequenceDictionary())
                .setOutputFile(config.OutputVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();

        addDepthHeaderLines(config, header);

        mWriter.writeHeader(header);

        mPendingRegions = new Semaphore(maxPendingRegions);
        mPendingVariants = Maps.newTreeMap();
        mNextWriteIndex = 0;
        mWriteLock = new Object();
    }

    public void acquireRegionPermit()
    {
        mPendingRegions.acquireUninterruptibly();
    }

    // returns a permit that was acquired but not used for a region (no task was left in the queue)
    public void releaseRegionPermit()
    {
        mPendingRegions.release();
    }

    // deposit a completed region's variants and flush now-contiguous regions to the VCF in order
    public void writeRegion(int regionIndex, final List<VariantContext> variants)
    {
        synchronized(mWriteLock)
        {
            mPendingVariants.put(regionIndex, variants);

            while(mPendingVariants.containsKey(mNextWriteIndex))
            {
                List<VariantContext> regionVariants = mPendingVariants.remove(mNextWriteIndex);
                regionVariants.forEach(x -> mWriter.add(x));
                ++mNextWriteIndex;
                mPendingRegions.release();
            }
        }
    }

    public void close()
    {
        if(!mPendingVariants.isEmpty())
            SV_LOGGER.warn("depth writer closed with {} pending regions unwritten",
                    mPendingVariants.size());

        mWriter.close();
    }

    private static void addDepthHeaderLines(final DepthConfig config, final VCFHeader header)
    {
        if(!header.hasFormatLine(ALLELE_FRACTION))
            header.addMetaDataLine(new VCFFormatHeaderLine(
                    ALLELE_FRACTION, 1, VCFHeaderLineType.Float, ALLELE_FRACTION_DESC));

        if(!header.hasFormatLine(REF_DEPTH))
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(
                    REF_DEPTH, 1, VCFHeaderLineType.Integer, REF_DEPTH_DESC));
            header.addMetaDataLine(new VCFInfoHeaderLine(
                    REF_DEPTH, 1, VCFHeaderLineType.Integer, REF_DEPTH_DESC));
        }
        else if(config.VcfTagPrefix != null)
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(
                    config.getVcfTag(REF_DEPTH), 1, VCFHeaderLineType.Integer, REF_DEPTH_DESC));
            header.addMetaDataLine(new VCFInfoHeaderLine(
                    config.getVcfTag(REF_DEPTH), 1, VCFHeaderLineType.Integer, REF_DEPTH_DESC));
        }

        if(!header.hasFormatLine(REF_DEPTH_PAIR))
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(
                    REF_DEPTH_PAIR, 1, VCFHeaderLineType.Integer, REF_DEPTH_PAIR_DESC));
            header.addMetaDataLine(new VCFInfoHeaderLine(
                    REF_DEPTH_PAIR, 1, VCFHeaderLineType.Integer, REF_DEPTH_PAIR_DESC));
        }
        else if(config.VcfTagPrefix != null)
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(
                    config.getVcfTag(REF_DEPTH_PAIR), 1,
                    VCFHeaderLineType.Integer, REF_DEPTH_PAIR_DESC));

            header.addMetaDataLine(new VCFInfoHeaderLine(
                    config.getVcfTag(REF_DEPTH_PAIR), 1,
                    VCFHeaderLineType.Integer, REF_DEPTH_PAIR_DESC));
        }
    }
}
