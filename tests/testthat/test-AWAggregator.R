test_that('Whether FragPipe aggregation workflow gives us the same output', {

    # From Workflow Ex.1
    data(sample.PSM.FP)
    regr <- loadPretrainedModel(useAvgCV=FALSE)
    samples <- c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9')
    groups <- samples
    df <- getAttributes(
        PSM=sample.PSM.FP,
        # TMT tag (229.1629) and carbamidomethylation (57.0214) are applied as fixed post-translational modifications (PTMs)
        fixedPTMs=c('229.1629', '57.0214'),
        colOfReporterIonInt=samples,
        groups=groups,
        setProgressBar=TRUE
    )
    aggregated_results <- aggregateByAttributes(
        PSM=df,
        colOfReporterIonInt=samples,
        ranger=regr,
        ratioCalc=FALSE
    )

    expect_equal(aggregated_results[1, 'Protein', drop = T], 'sp|A0AV96|RBM47_HUMAN')
    expect_equal(round(aggregated_results[2, 'Sample 2', drop = T], 4), 0.6601)
    expect_equal(round(aggregated_results[3, 'Sample 3', drop = T], 4), 1.0482)
})

test_that('Whether Proteome Discoverer aggregation workflow gives us the same output', {

    # From Workflow Ex.2
    data(sample.PSM.PD)
    data(sample.prot.PD)
    regr <- loadPretrainedModel(useAvgCV=FALSE)
    samples <- c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9')
    groups <- samples
    df <- convertPDFormat(
        PSM=sample.PSM.PD,
        protein=sample.prot.PD,
        colOfReporterIonInt=samples
    )
    df <- getAttributes(
        PSM=df,
        # TMT tag and carbamidomethylation are applied as static PTMs
        fixedPTMs=c('TMT6plex', 'Carbamidomethyl'),
        colOfReporterIonInt=samples,
        groups=groups,
        setProgressBar=TRUE
    )
    aggregated_results <- aggregateByAttributes(
        PSM=df,
        colOfReporterIonInt=samples,
        ranger=regr,
        ratioCalc=FALSE
    )

    expect_equal(aggregated_results[1, 'Protein', drop = T], 'A0AV96_Homo sapiens')
    expect_equal(round(aggregated_results[2, 'Sample 2', drop = T], 4), 0.6534)
    expect_equal(round(aggregated_results[3, 'Sample 3', drop = T], 4), 1.0495)
})

test_that('Whether getAttributes and getAvgScaledErrorOfLog2FC functions give us the same output', {

    # From Workflow Ex.3, step 2
    library(ExperimentHub)
    eh <- ExperimentHub()
    benchmarkSet3 <- eh[['EH9639']]
    samples <- c('H1+Y0_1', 'H1+Y1_1', 'H1+Y5_1', 'H1+Y10_1', 'H1+Y10_2', 'H1+Y5_2', 'H1+Y1_2', 'H1+Y0_2')
    groups <- c('H1+Y0', 'H1+Y1', 'H1+Y5', 'H1+Y10', 'H1+Y10', 'H1+Y5', 'H1+Y1', 'H1+Y0')
    PSM3 <- getAttributes(
        PSM=benchmarkSet3,
        fixedPTM=c('304.2071', '125.0476'),
        colOfReporterIonInt=samples,
        groups=groups,
        # The signals for yeast PSMs in group H1+Y0 is completely from noise, so they are not used for calculating Average CV
        groupsExcludedFromCV='H1+Y0'
    )
    PSM3 <- getAvgScaledErrorOfLog2FC(
        PSM=PSM3,
        colOfReporterIonInt=samples,
        groups=groups,
        expectedRelativeAbundance=list(`H1+Y0`=0, `H1+Y1`=1, `H1+Y5`=5, `H1+Y10`=10),
        speciesAtConstLevel='HUMAN'
    )

    expect_equal(round(PSM3[2, 'Scaled Intensity', drop = T], 4), 0.8016)
    expect_equal(round(PSM3[1, 'Distance Metric', drop = T], 4), 0.0793)
    expect_equal(round(PSM3[3, 'Average Scaled Error of log2FC', drop = T], 4), 0.6512)
})

