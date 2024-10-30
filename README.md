---
title: AWAggregator
author: Jiahua Tan, Gian L. Negri, Gregg B. Morin, David D. Y. Chen
output: 
  html_document: 
    keep_md: yes
    self_contained: no
---


### Introduction

The `AWAggregator` package implements an attribute-weighted aggregation algorithm which leverages peptide-spectrum match (PSM) attributes to provide a more accurate estimate of protein abundance compared to conventional aggregation methods. This algorithm employs pre-trained random forest models to predict the quantitative inaccuracy of PSMs based on their attributes. PSMs are then aggregated to the protein level using a weighted average, taking the predicted inaccuracy into account. Additionally, the package allows users to construct their own training sets that are more relevant to their specific experimental conditions if desired.

If you find this package helpful, please cite: 
XXXXXXXXXX

<br>


### Overview of Package Functions 

Functions available in the `AWAggregator` package:

* `get_dist_metric()`: Calculates the distance metric for PSMs.

* `get_attributes()`: Retrieves attributes required for training or test sets.

* `get_avg_scaled_error_of_log2FC()`: Calculates the Average Scaled Error of log2FC values required for training sets.

* `merge_training_sets()`: Extracts a similar number of PSMs from each input dataset and merges them into a single training set.

* `fit_model()`: Trains a random forest model to predict the level of quantitative inaccuracy of PSMs.

* `load_pretrained_model()`: Loads a pre-trained random forest model for predicting the level of quantitative inaccuracy of PSMs.

* `aggregate_by_attributes()`: Aggregates PSMs using a random forest model.

* `convert_PD_format()`: Converts output from Proteome Discoverer into the input format required by `AWAggregator`.

<br>

### Overview of Package Data

Data available in the `AWAggregator` package:

* `regr`, `child.nodeIDs`, `split.values.1`, `split.values.2`, `split.varIDs`: represent different parts of a pre-trained random forest model that incorporates the average coefficient of variation (CV) as a feature.

* `regr.no.CV`, `child.nodeIDs.no.CV`, `split.values.no.CV.1`, `split.values.no.CV.2`, `split.varIDs.no.CV`: represent different parts of a pre-trained random forest model that does not include the average CV as a feature.

* `benchmark.set.1`, `benchmark.set.2`, `benchmark.set.3`: represents PSMs in Benchmark Set 1 ~ 3 derived from the `psm.tsv` output files generated by FragPipe, which are used to train the random forest model. Columns unnecessary for the `AWAggregator` have been removed from the sample data.

* `sample.PSM.FP`: represents sample PSMs mapped to the proteins A0AV96, A0AVF1, A0AVT1, A0FGR8, and A0M8Q6, obtained from the `psm.tsv` output file generated by FragPipe. Columns unnecessary for the `AWAggregator` have been removed from the sample data.

* `sample.prot.PD`: represents sample proteins A0AV96, A0AVF1, A0AVT1, A0FGR8, and A0M8Q6, obtained from the TXT export of the proteins page in the Proteome Discoverer search results. Columns unnecessary for the `AWAggregator` have been removed from the sample data.

* `sample.PSM.PD`: represents sample PSMs mapped to the proteins A0AV96, A0AVF1, A0AVT1, A0FGR8, and A0M8Q6, obtained from the TXT export of the PSMs page in the Proteome Discoverer search results. Columns unnecessary for the `AWAggregator` have been removed from the sample data.

<br>


### Installation

The `AWAggregator` package can be installed using the `devtools` package.


```r
install.packages('devtools')
library(devtools)
install_github("Tan-Jiahua/AWAggregator")
```

<br>


### Workflow Examples

Load the `AWAggregator` package.


```r
library(AWAggregator)
```


<br>


#### Ex.1: Aggregate PSMs Derived from FragPipe Using the Pre-Trained Model.

In this example, we aggregate the reporter ion intensities of PSMs to the protein level. We use the sample dataset `sample.PSM.FP`, included in the `AWAggregator` package and derived from the `psm.tsv` output file generated by FragPipe. This dataset includes reporter ion intensities from nine samples, labeled from `Sample 1` to `Sample 9`, without replicates. The PSMs are mapped to the following proteins: A0AV96, A0AVF1, A0AVT1, A0FGR8, and A0M8Q6, with unnecessary columns removed for clarity.

This example demonstrates the basic functionality of the `AWAggregator` package using the default pre-trained model.


```r
# Load the pre-trained random forest model that does not include the average CV as a feature, which indicates the average CV in percentage for processed PSM reporter ion intensities across different replicate groups. It is recommended to load the pre-trained model with average CV when replciates are available; otherwise, use the model without the average CV
regr = load_pretrained_model(use_average_CV = F)
df = get_attributes(
  PSM = sample.PSM.FP,
  # TMT tag (229.1629) and carbamidomethylation (57.0214) are applied as fixed post-translational modifications (PTMs)
  fixed_PTMs = c('229.1629', '57.0214'),
  col_of_reporter_ion_int = c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9'),
  groups = c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9'),
  set_progress_bar = T
)
aggregated_results = aggregate_by_attributes(
  PSM = df,
  col_of_reporter_ion_int = c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9'),
  ranger = regr
)
```


<br>

The output dataframe will provide estimates of protein abundance.


```
Protein               Sample 1   Sample 2   Sample 3   Sample 4   ...
sp|A0AV96|RBM47_HUMAN 0.11216715 0.12205447 0.09577089 0.11595996 ...
sp|A0AVF1|IFT56_HUMAN 0.08646259 0.08586356 0.08711255 0.10386396 ...
sp|A0AVT1|UBA6_HUMAN  0.12710276 0.12570251 0.11212039 0.11669525 ...
sp|A0FGR8|ESYT2_HUMAN 0.10666391 0.09749581 0.12095627 0.09117375 ...
sp|A0M8Q6|IGLC7_HUMAN 0.05410301 0.08225112 0.09616674 0.08964740 ...
```


<br>


#### Ex.2: Aggregate PSMs Derived from Proteome Discoverer Using the Pre-Trained Model.

In this example, we convert the search result from Proteome Discoverer to the format required by `AWAggregator` and aggregate the reporter ion intensities of PSMs to the protein level. We use the sample dataset `sample.PSM.PD`, alongside its corresponding protein table `sample.prot.PD`, both included in the `AWAggregator` package. These files are derived from the TXT exports of the proteins and PSMs pages in the search results from Proteome Discoverer. This dataset includes reporter ion intensities from nine samples, labeled from `Sample 1` to `Sample 9`, without replicates. The PSM and protein tables contains following proteins: A0AV96, A0AVF1, A0AVT1, A0FGR8, and A0M8Q6, with unnecessary columns removed for clarity.


```r
# Load the pre-trained random forest model that does not include the average CV as a feature, which indicates the average CV in percentage for processed PSM reporter ion intensities across different replicate groups. It is recommended to load the pre-trained model with average CV when replciates are available; otherwise, use the model without the average CV
regr = load_pretrained_model(use_average_CV = F)
df = convert_PD_format(
  PSM = sample.PSM.PD,
  protein = sample.prot.PD,
  col_of_reporter_ion_int = c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9')
)
df = get_attributes(
  PSM = df,
  # TMT tag and carbamidomethylation are applied as static PTMs
  fixed_PTMs = c('TMT6plex', 'Carbamidomethyl'),
  col_of_reporter_ion_int = c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9'),
  groups = c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9'),
  set_progress_bar = T
)
aggregated_results = aggregate_by_attributes(
  PSM = df,
  col_of_reporter_ion_int = c('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9'),
  ranger = regr
)
```


<br>

The output dataframe will provide estimates of protein abundance.

```
Protein               Sample 1   Sample 2   Sample 3   Sample 4   ...
1 A0AV96_Homo sapiens 0.11848943 0.12003883 0.08952637 0.11850774 ...
2 A0AVF1_Homo sapiens 0.08119006 0.08048804 0.08771501 0.09536281 ...
3 A0AVT1_Homo sapiens 0.12526313 0.12122091 0.10922527 0.11575040 ...
4 A0FGR8_Homo sapiens 0.10660835 0.09256367 0.12074519 0.08639497 ...
5 A0M8Q6_Homo sapiens 0.04614932 0.06161329 0.09481033 0.07929251 ...
```


<br>


#### Ex.3: Build a Merged Training Set and Retrain the Model.

Retraining the AWA model using additional spike-in datasets can improve the number of quantified PSMs in the merged training set, and hence the robustness of the correlation. In addition, retraining using experiment-specific in-house spike-in datasets could also provide potential benefits for the machine learning model by better representing the employed hardware and acquisition modes.

In this example, we create a training set by merging three benchmark spike-in datasets (`benchmark.set.1`, `benchmark.set.2`, and `benchmark.set.3`), all included in the `AWAggregator` package and derived from the `psm.tsv` output files generated by FragPipe. This combined training set is then used to train a random forest model.


<br>

Firstly, we calculate the attributes and the values of Average Scaled Error of log<sub>2</sub>FC in `benchmark.set.1`.

```r
samples = c('H1+E1_1', 'H1+E1_2', 'H1+E1_3', 'H1+E2_1', 'H1+E2_2', 'H1+E2_3', 'H1+E2_4', 'H1+E6_1', 'H1+E6_2', 'H1+E6_3')
groups = c('H1+E1', 'H1+E1', 'H1+E1', 'H1+E2', 'H1+E2', 'H1+E2', 'H1+E2', 'H1+E6', 'H1+E6', 'H1+E6')
PSM1 = get_attributes(
  PSM = benchmark.set.1,
  # TMT tag (229.1629) and carbamidomethylation (57.0214) are applied as fixed PTMs
  fixed_PTM = c('229.1629', '57.0214'),
  col_of_reporter_ion_int = samples,
  groups = groups
)
PSM1 = get_avg_scaled_error_of_log2FC(
  PSM = PSM1,
  col_of_reporter_ion_int = samples,
  groups = groups,
  # The actual protein fold change may be deviated from the intended values after TMT labelling as the original work indicates when H1+Y6 is involved, and therefore, H1+Y6 is not used in the calculation of Average of Scaled Error of log2FC
  expected_relative_abundance = list(`H1+E1` = 1, `H1+E2` = 2, `H1+E6` = NA),
  species_at_const_level = 'HUMAN'
)
```

<br>

Secondly, we calculate the attributes and the values of Average Scaled Error of log<sub>2</sub>FC in `benchmark.set.2`. `benchmark.set.2` consists of three separate mass spectrometry runs, each of which is processed individually.

```r
library(dplyr)
samples = c('H1+Y1_1', 'H1+Y4_1', 'H1+Y10_1', 'H1+Y1_2', 'H1+Y4_2', 'H1+Y10_2', 'H1+Y1_3', 'H1+Y4_3', 'H1+Y10_3')
groups = c('H1+Y1', 'H1+Y4', 'H1+Y10', 'H1+Y1', 'H1+Y4', 'H1+Y10', 'H1+Y1', 'H1+Y4', 'H1+Y10')
PSM2 = bind_rows(lapply(unique(benchmark.set.2$Replicate), FUN = function(X){
  df = get_attributes(
    PSM = benchmark.set.2[benchmark.set.2$Replicate == X, ],
    fixed_PTM = c('229.1629', '57.0214'),
    col_of_reporter_ion_int = samples,
    groups = groups,
    set_progress_bar = F
  )
  df = get_avg_scaled_error_of_log2FC(
    PSM = df,
    col_of_reporter_ion_int = samples,
    groups = groups,
    expected_relative_abundance = list(`H1+Y1` = 1, `H1+Y4` = 4, `H1+Y10` = 10),
    species_at_const_level = 'HUMAN'
  )
  return(df)
}))
```

<br>

Thirdly, we calculate the attributes and the values of Average Scaled Error of log<sub>2</sub>FC in `benchmark.set.3`.

```r
samples = c('H1+Y0_1', 'H1+Y1_1', 'H1+Y5_1', 'H1+Y10_1', 'H1+Y10_2', 'H1+Y5_2', 'H1+Y1_2', 'H1+Y0_2')
groups = c('H1+Y0', 'H1+Y1', 'H1+Y5', 'H1+Y10', 'H1+Y10', 'H1+Y5', 'H1+Y1', 'H1+Y0')
PSM3 = get_attributes(
  PSM = benchmark.set.3,
  # TMTpro tag (304.2071) and N-ethylmaleimide (125.0476) are applied as fixed PTMs
  fixed_PTM = c('304.2071', '125.0476'),
  col_of_reporter_ion_int = samples,
  groups = groups,
  # The signals for yeast PSMs in group H1+Y0 is completely from noise, so they are not used for calculating Average CV
  groups_excluded_from_CV = 'H1+Y0'
)
PSM3 = get_avg_scaled_error_of_log2FC(
  PSM = PSM3,
  col_of_reporter_ion_int = samples,
  groups = groups,
  expected_relative_abundance = list(`H1+Y0` = 0, `H1+Y1` = 1, `H1+Y5` = 5, `H1+Y10` = 10),
  species_at_const_level = 'HUMAN'
)
```

<br>

Next, we merge a new training set from these three datasets. The minimum number of PSMs to extract from each dataset is determined by the number of PSMs in the smallest set. Complete sets of PSMs mapped to the selected proteins are extracted, resulting in final PSM counts from each set that are equal to or slightly larger than the preset values.

```r
PSM = merge_training_sets(
  PSM_list = list(
    `Benchmark Set 1` = PSM1,
    `Benchmark Set 2` = PSM2,
    `Benchmark Set 3` = PSM3
  ),
  num_PSMs = min(nrow(PSM1), nrow(PSM2), nrow(PSM3)),
  seed = 1000
)
```

<br>

Train a new random forest model using Average CV as an attribute.

```r
regr = fit_model(PSM, use_average_CV = T, seed = 3979)
```
