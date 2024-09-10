
# GSAuto - Methylation IDAT File Processing

## Overview

**GSAuto** is an R function designed to process IDAT files from Illumina methylation arrays (27k, 450k, 850k, 935k, MSA). It extracts methylation intensity data, calculates beta values, p-values, and other quality metrics, and organizes the results into matrices.

## Features

- Supports multiple platforms: Illumina 27k, 450k, 850k(EPIC v1.0), 935k(EPIC v2.0), and MSA.
- Processes `.idat` files and calculates:
  - **Beta values**: The methylation level.
  - **P-values**: Derived from comparison with control probes.
  - **Detection rates**: Reflecting the quality of the sample based on p-values.
  - **Intensity values**: Summed signal intensities for both methylated and unmethylated probes.

## Prerequisites

Before using the function, ensure that you have the following R packages installed:

```r
install.packages("illuminaio")
install.packages("dplyr")
install.packages("data.table")
```

Additionally, make sure you have access to the manifest and control files for the respective platform stored in the appropriate directories.

## Installation

Download the R script directly from this repository. After downloading, load it into your R session:

```r
source("GSAuto.R")
```

## Function Usage

```r
result_list <- GSAuto(idat_dir = "/path/to/idat/files", platform = "450k")
```

### Parameters:
- **idat_dir**: The directory where your `.idat` files are stored.
- **platform**: The platform type (e.g., "27k", "450k", "850k", "935k", "MSA").

### Output:
The function returns a list containing the following elements:
- `detect_re`: A data frame of detection rates and quality control results for each sample.
- `beta_matrix`: A matrix of beta values for each sample.
- `pval_matrix`: A matrix of p-values for each sample.
- `M_matrix`: Methylated intensities for each sample.
- `U_matrix`: Unmethylated intensities for each sample.
- `Intensity_matrix`: Total intensities for each sample.

## Example

```r
# Example usage
result_list <- GSAuto(idat_dir = "/path/to/idat", platform = "450k")

# Access the beta values matrix
beta_matrix <- result_list$beta_matrix

# Access the p-value matrix
pval_matrix <- result_list$pval_matrix
```

## File Structure

The directory should contain the following files based on the platform you are using:

```text
/path/to/data/gsauto_manifest/
  ├── 27k/
  │   ├── negControl.rds
  │   ├── manifest.rds
  ├── 450k/
  │   ├── negControl.rds
  │   ├── manifest.rds
  ├── 850k/
  │   ├── negControl.rds
  │   ├── manifest.rds
  ├── 935k/
  │   ├── negControl.rds
  │   ├── manifest.rds
  ├── MSA/
      ├── negControl.rds
      ├── manifest.rds
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
