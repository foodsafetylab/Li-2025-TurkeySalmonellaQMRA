# Li-2025-TurkeySalmonellaQMRA

## Overview
Data and code used for paper titled "Risk Assessment Predicts Final Product Standards Most Specifically Targeting Products with High-Level of Highly-Prevalent High-Virulence Serotypes Most Efficiently Reduce Risk of *Salmonella* from Comminuted Turkey" can be found here. 

In 2024, USDA-FSIS published risk assessments for Salmonella in raw poultry products. Their risk assessments supported final product standards targeting high-level Salmonella and highly prevalent high-virulence serotypes. To assess the impact of more specific standards, we adapted our previous risk assessment for chicken parts to comminuted turkey, incorporating critical assumptions from FSIS risk assessments. We simulated the attributable fraction of illnesses from products contaminated over three level thresholds (0.0031 CFU/g, 1 CFU/g, and 10 CFU/g) and/or containing a serotype in three lists (“Top 3 most prevalent high-virulence serotypes”, “All high-virulence serotypes”, and “High-virulence proportion of each serotype”). Results showed that 87% of illnesses were attributed to the 0.56% of products with Salmonella exceeding 10 CFU/g. Under more specific criteria of level “AND” serotype, 60% of illnesses were attributed to the 0.14% of products contaminated with Salmonella exceeding 10 CFU/g and one of the three most prevalent high-virulence serotypes. Further, applying genomic-based clustering information, 75% of illnesses were attributed to slightly more products (0.19%) containing Salmonella exceeding 10 CFU/g and high-virulence proportion of each serotype. Under the less specific standard, however, 99% of illnesses were attributed to the 5.7% of products containing Salmonella exceeding 10 CFU/g “OR” one of the three most prevalent high-virulence serotypes. Our study demonstrated that most salmonellosis risk is concentrated in comminuted turkey products with high levels of high-virulence contaminations. Importantly, more specifically targeting those products could efficiently reduce public health risk while minimizing products implicated.
## Usage
### Setup
- **Raw data** was extracted from USDA-FSIS HACCP verification from [FSIS Laboratory Sampling Data Webpage](https://www.fsis.usda.gov/news-events/publications/raw-poultry-sampling) accessed September 16, 2024. Raw data can be found [here](/Raw%20data). 
- ***Salmonella* level fitting** and **Serotype distribution** for both 2016-2021 and 2023-2024 datasets can be found [here](/Level%20fitting%20and%20Serotype%20distribution).
- **Risk assessment models** for the QMRA can be found [here](/Risk%20model). 

### Running
**Risk assessment models** can be run by opening the Risk Model.Rproj file, followed by "Turkey model_baseline.Rmd" for scenarios using the majority clustering method and "Turkey model_proportion.Rmd" for scenarios using the proportion clustering method. Please follow the annotations and run the required chunks of code. 

## Authors
You can view the list of authors in the [AUTHORS](/AUTHORS) file.

## Contact
Corresponding author: Matthew J. Stasiewicz<br>
103 Agricultural Bioprocess Lab<br>
1302 W. Pennsylvania<br>
Urbana, IL, 1361801<br>
USA<br>
+1-217-265-0963<br>
[mstasie@illinois.edu](mailto:mstasie@illinois.edu)

## Citation
Include citation here.

## License
This project's code is licensed under the GNU General Public License v3.0 and dataset is licensed the Creative Commons Attribution Share Alike 4.0 International license. Please see the [LICENSE.code](/LICENSE.code) and [LICENSE.dataset](/LICENSE.dataset) files for details.

## Funding
This study was an extension project of #BRF-015 Risk Assessment Comparing Alternative Approaches to Regulating Salmonella in Poultry by Public Health Impact Factors, funded by US Poultry and Egg Association. This work was supported by the Jonathan Baldwin Turner Fellowship, awarded to Yiyi Li, from the College of Agricultural, Consumer and Environmental Sciences at the University of Illinois Urbana-Champaign. 
