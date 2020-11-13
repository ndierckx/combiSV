# combiSV

Combine structural variation outputs from long sequencing reads into a superior call set

### Getting help

Any issues/requests/problems/comments that are not yet addressed on this page can be posted on [Github issues](https://github.com/ndierckx/Sim-it/issues) and I will try to reply the same day.

Or you can contact me directly through the following email address:

nicolasdierckxsens at hotmail dot com 

### Prerequisites

Perl<br> 

### Instructions

Usage:

<code>perl combiSV1.0.pl -pbsv <pbsv_output.vcf> -sniffles <sniffles_output.vcf> -nanovar <nanovar_output.vcf> -svim <svim_output.vcf> -nanosv <nanosv_output.vcf> -c <minimum_variance_allele_coverage> -o <output_name></code>
 

### Output

#### 1. output_name.vcf: 
This is the combined standard vcf output 

#### 2. simplified_output_name.vcf: 
This is a simplified vcf output that can be used as input for Sim-it
