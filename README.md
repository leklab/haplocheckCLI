Refactoring of code https://github.com/genepi/haplocheck such that it takes vcf file containing directory from command line outside the cloudgene framework

**Compile**  
jar cvfe haplocheckCLI.jar haplocheck_contam *

**Usage**  
java -jar haplocheckCLI.jar directory_with_vcf

**TO DO**  
Switch back from VL to AF  
Get rid of those annoying double quotes in output  
Other cool stuff
