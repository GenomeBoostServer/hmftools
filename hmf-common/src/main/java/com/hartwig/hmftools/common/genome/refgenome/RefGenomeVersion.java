package com.hartwig.hmftools.common.genome.refgenome;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum RefGenomeVersion
{
    V37("37", true, false),
    V38("38", false, true),
    // hg19: GRCh37 coordinates (is37) but chr-prefixed chromosome naming (like V38)
    V19("19", true, true);

    @NotNull
    private final String mIdentifier;
    private final boolean mIs37;
    private final boolean mChrPrefixed;

    // config option
    public static final String REF_GENOME_VERSION = "ref_genome_version";
    public static final String REF_GENOME_VERSION_CFG_DESC =
            "Ref genome version: 37, 38 or 19 (hg19 = GRCh37 coords with chr prefix)";

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeVersion.class);
    private static final String GZIP_EXTENSION = ".gz";

    @NotNull
    public static RefGenomeVersion from(@NotNull final String version)
    {
        if(version.equalsIgnoreCase(V37.toString()) || version.equalsIgnoreCase("37")
                || version.equalsIgnoreCase("HG37"))
        {
            return V37;
        }
        else if(version.equalsIgnoreCase(V38.toString()) || version.equalsIgnoreCase("38")
                || version.equalsIgnoreCase("HG38"))
        {
            return V38;
        }
        else if(version.equalsIgnoreCase(V19.toString()) || version.equalsIgnoreCase("19")
                || version.equalsIgnoreCase("HG19"))
        {
            return V19;
        }

        throw new IllegalArgumentException("Cannot resolve ref genome version: " + version);
    }

    public static RefGenomeVersion from(final CommandLine cmd) { return from(cmd, V37); }

    public static RefGenomeVersion from(final CommandLine cmd, final RefGenomeVersion defaultVersion)
    {
        return RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION, String.valueOf(defaultVersion)));
    }

    RefGenomeVersion(
            @NotNull final String identifier, final boolean is37, final boolean chrPrefixed)
    {
        mIdentifier = identifier;
        mIs37 = is37;
        mChrPrefixed = chrPrefixed;
    }

    public boolean isChrPrefixed() { return mChrPrefixed; }

    public boolean is37() { return mIs37; }
    public boolean is38 () { return !mIs37; }

    public String identifier() { return mIdentifier; }

    @NotNull
    public String versionedChromosome(@NotNull String chromosome)
    {
        // chr-prefix behaviour is driven by naming convention (V38 and V19 are chr-prefixed),
        // not by the coordinate version
        return mChrPrefixed ?
                RefGenomeFunctions.enforceChrPrefix(chromosome) :
                RefGenomeFunctions.stripChrPrefix(chromosome);
    }

    public String addVersionToFilePath(final String filePath)
    {
        String modifiedFilePath = filePath;
        if(filePath.endsWith(GZIP_EXTENSION))
        {
            modifiedFilePath = filePath.substring(0, filePath.indexOf(GZIP_EXTENSION));
        }

        if(!modifiedFilePath.contains("."))
        {
            throw new IllegalStateException("Cannot include ref genome version in file path that has no proper extension: " + filePath);
        }

        int extensionStart = modifiedFilePath.lastIndexOf(".");
        String versionedFilePath =
                modifiedFilePath.substring(0, extensionStart) + "." + this.mIdentifier + modifiedFilePath.substring(extensionStart);

        if(filePath.endsWith(GZIP_EXTENSION))
        {
            versionedFilePath = versionedFilePath + GZIP_EXTENSION;
        }

        return versionedFilePath;
    }
}
