# MCuptake

Package for estimating male circumcision uptake rates

# Introduction

Male circumcision reduces the risk of HIV acquisition in men. HIV impact models like the [Spectrum Goals](https://avenirhealth.org/software-spectrum.php) model use estimates of male circumcision prevalence to estimate and project HIV incidence. This package was designed to help estimate male circumcision prevalence trends by five-year age group using data from nationally-representative household surveys conducted in countries prioritized for scale-up of voluntary medical male circumcision (VMMC).

# How to use the MCuptake package

In R, you can install the package from GitHub using `devtools`

```         
devtools::install_git('rlglaubius/MCuptake')
```

Once the package is installed, you can run the `tester.R` script in this repository to run an example fit. This script is configured to estimate male circumcision prevalence in Zimbawbwe out-of-the-box. To estimate trends for a different country, change the three-character country code in the script from "ZWE" to one of the other supported countries (see table in the `Data` section below).

# Data

The MCuptake package uses population size estimates by year and age from the [United Nations Population Division 2022 revision of the World Population Prospects](https://population.un.org/wpp/).

Household survey data were obtained from Demographic Health Surveys (DHS), AIDS Indicator Surveys (AIS), and Population-Based HIV Impact Assessments (PHIA). Data were obtained from published survey reports or from [StatCompiler](https://www.statcompiler.com). After loading the `MCuptake` package, the survey data is stored in the variable `mc_svy_data`.

| Country                     | Country Code |
|-----------------------------|--------------|
| Eswatini                    | SWZ          |
| Ethiopia                    | ETH          |
| Kenya                       | KEN          |
| Lesotho                     | LSO          |
| Malawi                      | MWI          |
| Mozambique                  | MOZ          |
| Namibia                     | NAM          |
| Rwanda                      | RWA          |
| South Africa                | ZAF          |
| Uganda                      | UGA          |
| United Republic of Tanzania | TZA          |
| Zambia                      | ZMB          |
| Zimbabwe                    | ZWE          |

Note that two priority countries for VMMC scale-up, Botswana and South Sudan, are not included. Since Botswana's 2021 BAIS V survey included male circumcision prevalence, we anticipate adding Botswana in a future update.
