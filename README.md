# Sv6_StarlingEcogeographical

This repository contains code related to the following publication:

Stuart KC, Sherwin WB, Cardilini APA, Rollins LA **2022**. Genetics and Plasticity Are Responsible for Ecogeographical Patterns in a Recent Invasion. *Frontiers in Genetics*, 13: 824424, [10.3389/fgene.2022.824424](10.3389/fgene.2022.824424)

In the raw scripts folder you will find:
* 2022_02_12_notebook_16949.pdf: GCTA analysis
* 2022_02_12_notebook_16923.pdf: Morphology-associated alleles correlation analysis (also covered below in the vignette)
* 2022_02_12_notebook_16913.pdf: GradientForest analysis

<h3>Morphology-associated alleles correlation analysis</h3>

The main aim of this project was to explore ways that you can combine trio data sets of genetic, phenotype (morphology), and spatial (enviornmental) data. One of the ways in which I did this was to decide to focus on morphology-associated SNPs from my data set and comparing these to the environmental data I had from sampling site coordinates, to obtain portential climate drivers of the genetics underlying morphology patterns. This was the first time I had oppertunity to do such analysis (morphology data paired with genetic data is surpringly hard to come by in animal invasion genomics), so it could probably use some refining in future project. But for now, I walk through the code I have used below.

The general analysis pipeline is:
* Format files for BayPass analysis in plink
* Find morphology-associated SNPs in BayPass
* Produce reduced VCF files that contain only these SNPs, and summarise with prinicple component analysis
* Regress these against climate variables, calculate contribution of each climate var to the model

This analysis was essentially run separately for each morphological variable: mass, tarsus (length), head (width), beak (surface area), wing (length), spleen (mass), and heart (mass).

## Software

Below is a list of the modules used (and versions).

<pre class="r"><code>module load stacks/2.2
module add vcftools/0.1.16
module load plink/1.90b6.7
module load baypass
module load R/3.6.3
module load bedtools/2.27.1
</code></pre>

## Plink file processing (for BayPass analysis)

Please refer to the [BayPass Manual](http://www1.montpellier.inra.fr/CBGP/software/baypass/files/BayPass_manual_2.3.pdf); here we are aiming to produce a 'genotyping data file'. We start with a file in VCF format, and use Stacks populations to convert this to plink file format.

<pre class="r"><code>VCF=australian_starlings.vcf
populations -V $VCF -O ./ -M australian_starlings_popmap.txt -t 2 --plink
</code></pre>

Now we will use plink to summarise genotype allele data for each population. As an extra note, each morphological variable will be run separately in BayPass. Because of this, and because I had incomplete morphology data set (some individuals were missing morphology data for one measure, however had all other measures) I had to kick out some individuals for each separate morphology analysis (these are captured by 'GCTA_${PHENO}_subset.txt'). Reterospectively, I could have left these in as BayPass works on population averages, but I still think this is the more proper way of doing it. 

<pre class="r"><code>for PHENO in mass tarsus head beak wing spleen heart
do
plink --file australian_starlings_popmap.plink  --keep ../GCTA212_7pheno/GCTA_${PHENO}_subset.txt --allow-extra-chr --freq counts --family --out plink_files/plink_${PHENO}
done
</code></pre>

We're not done formatting yet. Here, we manipulate the plink output so that it contains the necessary data for a BayPass 'genotyping data file'. The genotyping data file is simply organized as a matrix with nsnp rows and 2 ∗ npop columns, which we can extract directly from the plink output. Please note, that the below one liner should work for your data, except you need to change the **28** in the sed command to be your own *population/sampling site/family* count multiplied by 2. I had 14 sample sites, so 14*2=28, therefore I add 28 here. 

A super brief explination of the below one liner is: tail to remove the plink file header | produce an allele count for major and minor allele | retain only the major and minor allele counts | clean up some of the characters | create a new line every 28 fields. 

<pre class="r"><code>for PHENO in mass tarsus head beak wing spleen heart
do
tail -n +2 plink_files/plink_${PHENO}.frq.strat | awk '{ $9 = $8 - $7 } 1' | awk '{print $7,$9}' | tr "\n" " " | sed 's/ /\n/28; P; D'> genofile_baypass_${PHENO}.txt
done
</code></pre>

The genotyping data file is now ready for BayPass.

## Running BayPass

From now on most of these were done as array scripts, so assume the code has been run separately for each morphological variable, with "PHENO" being substituted out for an appropriate variable/file name.

Run the base BayPass model. 

<pre class="r"><code>g_baypass -npop 14 -gfile ../genofile_baypass_${PHENO}.txt -outprefix Pheno_${PHENO}_anacore -nthreads 16
</code></pre>

Now we have the omegafile output we can run the association models. The 'baypass_pheno_${PHENO}.txt' just contains the average morphology measure for each population. 

Please note that I ran both the anacovaux and anacovis models - I really should have just picked one, and I think most people pick anacovaux. Sometimes running a bunch of different things within a program can help your learning when you are just starting out in genomics. My thresholds were pretty high for the SNPs I chose to retain after this analysis, so I don't think it was super important in terms of the result, but this is definitely a trap I fell into early on. There are just so many options and when you are new working with genomics data they can get a bit overwhelming and maing decisions feels hard/arbitrary. 

<pre class="r"><code>g_baypass -npop 14 -gfile ../genofile_baypass_${PHENO}.txt -efile ../baypass_pheno_${PHENO}.txt -scalecov -auxmodel -nthreads 16 -omegafile Pheno_${PHENO}_anacore_mat_omega.out -outprefix Pheno_${PHENO}_anacovaux_scaled

g_baypass -npop 14 -gfile ../genofile_baypass_${PHENO}.txt  -efile ../baypass_pheno_${PHENO}.txt -scalecov -nthreads 16 -omegafile Pheno_${PHENO}_anacore_mat_omega.out -outprefix Pheno_${PHENO}_anacovis_scaled
</code></pre>

Now we can plot the association with this morphology variable along the length of the genome. For this we jump out of the unix environment an into R.

Here is the code for the anacovis plot for the morphology measure mass.

<pre class="r"><code>R
setwd("/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv6_Morphology/analysis/baypass_phenofilter212_7pheno/Pheno_mass")
covis.snp.res.mass=read.table("Pheno_mass_anacovis_scaled_summary_betai_reg.out",h=T)
graphics.off()
pdf("Pheno_mass_covis_scaled.pdf")
layout(matrix(1:3,3,1))
plot(covis.snp.res.mass$BF.dB.,xlab="SNP",ylab="BFis (in dB)")
abline(h=20, col="red")
plot(covis.snp.res.mass$eBPis,xlab="SNP",ylab="eBPis")
plot(covis.snp.res.mass$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()
</code></pre>

![screenshot](/Sv6_vignette/Sv6_covis.png)

And for the  covaux plot.

<pre class="r"><code>R
setwd("/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv6_Morphology/analysis/baypass_phenofilter212_7pheno/Pheno_mass")
covaux.snp.res.mass=read.table("Pheno_mass_anacovaux_scaled_summary_betai.out",h=T)
covaux.snp.xtx.mass=read.table("Pheno_mass_anacovaux_scaled_summary_pi_xtx.out",h=T)$M_XtX
graphics.off()
pdf("Pheno_mass_covaux_scaled.pdf")
layout(matrix(1:3,3,1))
plot(covaux.snp.res.mass$BF.dB.,xlab="Mass",ylab="BFmc (in dB)")
abline(h=20, col="red")
plot(covaux.snp.res.mass$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"))
plot(covaux.snp.xtx.mass, xlab="SNP",ylab="XtX corrected for SMS")
dev.off()
</code></pre>

![screenshot](/Sv6_vignette/Sv6_covaux.png)

For the above two plots, I chose a BF of 20. This decision is a little arbitrary (all statistical thresholds are...), but I would recommend reading further about Jeffrey’s rule and interpreting BF outputs when you are setting your own threshold. You can even set a few and catagorise SNPs within [threshold groups](https://rdrr.io/cran/effectsize/man/interpret_bf.html).

## Filtering your SNPs for those above the BF threshold

Baypass numbers the SNP output, so here I add line numbers to the SNP list file so that I can match them up to the Baypass line numbers. NOTE: while we did remove some individuals that had missing data from the variant files, this did not impact the total SNP count per phenotype-specific VCF. Hence we can use the samme numbered reference VCF for each phenotype Baypass output.

Brief summary of this one liner: grab only non-header lines from VCF | keep only first 3 columns (the snp chrom, pos, and name) | remove header | print line number at the end of the row.

<pre class="r"><code>VCF=australian_starlings.vcf
grep -v "^##" {VCF} | cut -f1-3 | tail -n +2 | awk '{print $0,NR}' > australian_starlings_snplist_numbered.txt
</code></pre>

Now we can use the above file to work out which SNPs Baypass reported as morphology-associated (i.e. outliers). The below code pools the outliers from the anacovaux and the anacovis model, and retains a list of which SNPs these correspond to.

<pre class="r"><code>cd Pheno_${PHENO}/
cat Pheno_${PHENO}_anacovaux_scaled_summary_betai.out | awk '$6>20' > outliers_${PHENO}_anacovaux_BF20.txt
cat Pheno_${PHENO}_anacovis_scaled_summary_betai_reg.out | awk '$5>20' > outliers_${PHENO}_anacovis_BF20.txt
awk 'FNR==NR{a[$2];next} (($4) in a)' outliers_${PHENO}_anacovaux_BF20.txt australian_starlings_snplist_numbered.txt > outliers_${PHENO}_anacovaux_BF20_SNPlist.txt
awk 'FNR==NR{a[$2];next} (($4) in a)' outliers_${PHENO}_anacovis_BF20.txt australian_starlings_snplist_numbered.txt > outliers_${PHENO}_anacovis_BF20_SNPlist.txt
sort outliers_${PHENO}_anacovis_BF20_SNPlist.txt outliers_${PHENO}_anacovaux_BF20_SNPlist.txt | uniq > outliers_${PHENO}_combine_BF20_SNPlist.txt
</code></pre>

## Conducting PCA on morphology-associated snps

How we have a list of morphology-associated SNPs, we can return to the original VCF, filter this for just the morphology-associated SNPs, and conduct PCA. The last line of this 

<pre class="r"><code>cut -f3 ../Pheno_${PHENO}/outliers_${PHENO}_combine_BF20_SNPlist.txt > ${PHENO}_snps.txt
VCF=australian_starlings.vcf
vcftools --vcf ${VCF} --snps ${PHENO}_snps.txt --recode --out ${PHENO}_snps
vcftools --vcf ${PHENO}_snps.recode.vcf --plink --out plink_files/${PHENO}_snps.plink
plink --file plink_files/${PHENO}_snps.plink --pca --out plink_files/${PHENO} --make-rel
</code></pre>

## Regression with environmental variables

Next, we pull the data into R to conduct the regressions. First load up the packages we need. We also load in out meta data and filter for just the column information we need. You metadata file will of course look different from mine, but essentially you just need the individual ID and environmental/climate/spatial variables you are interested in.

<pre class="r"><code>library(dplyr)
library(psych)
library(data.table)
library(relaimpo)
Starling <- read.csv("australian_starlings_metadata.csv",stringsAsFactors=TRUE,sep=",")
Starling_pheno_212_env <- filter(Starling, PhenoDataTarsus_212 == "YES") %>% select(c(location, INDV, bio01, bio03, bio04, bio08, bio09, bio12, bio14, bio15, bio18, bio19, elev, mnNDVI))
str(Starling_pheno_212_env)
</code></pre>

I truely apologise for the next piece of code. I was trying to apply unix loop logic to R, and it got messy. I am sure the same thing can be achieved with a much neater R function, however at the time I was blisfully unaware that this was an option. I'll run through an example for one phenotype: mass.

NOTE: If I were to do this, I think the linear model should include some kind of spatial random variable to account for distance based environmental differences. If you yourself are trying to do something similar, feel free to reach out and discuss model logic/choice. 

Below, the genetic PCA axis eigenvalues serving as the response variable, and a selection of environmental variables serving as the predictors (variables selected using Gradient Forest analysis). Originally I planned to use all 20 PCA axes (hence the loop below for 1:20), and the overall model R2 values (amount of genetic PCA axis variance explained by environmental predictors) were weighted by the percentage variance each of the genetic axis captured of the overall variance across all genetic PCA axis. However I realised early on that this cut off was quite arbitrary - the first 20 axes accounted for quite different proportions of the overall genetic variation for each genetic data subset. The final 2 lines of code keep just the axis that explain above a threshold of variance. So instead I chose to retain only the axis above the plateau, as the assumption is that keeping only this first chunk of variance will allow me to assess if the major patterns of variation in the data set that appeared to be correlated with the environmental variables. These turned out to be somewhere between 1-3 axes for each morphology variable, but I left this code in there in case I needed to refer back to it in the future (and because it was useful for the morphology data that required more than 1 axis).

<pre class="r"><code>pca.mass.import <- read.table("plink_files/mass.eigenvec", sep=" ", header=F)
lmer.all.mass <- NA
lmer.all.mass.variance <- NA

for (i in 1:20) {
         j <- i + 2    
pca.mass <- pca.mass.import[,c(2,j)]
     names(pca.mass)[1] <- "INDV"
names(pca.mass)[2] <- "PCA"
modeldata.mass <- merge(Starling_pheno_212_env, pca.mass, by.x = "INDV")
fit.mass <- lm( PCA ~ bio01 + bio03 + bio04 + bio08 + bio09 + bio12 + bio14 + bio15 + bio18 + bio19 + elev + mnNDVI,data=modeldata.mass)
     imp.mass <- calc.relimp(fit.mass,type=c("lmg"), rela=TRUE)
lmer.all.mass <- data.frame(lmer.all.mass,imp.mass@lmg)
         z <- i + 1
names(lmer.all.mass)[z] <- paste0("PCA", i)
lmer.all.mass.variance <- data.frame(lmer.all.mass.variance,imp.mass@R2)
         z <- i + 1
names(lmer.all.mass.variance)[z] <- paste0("PCA", i)
}

pca.mass.vals <- read.table("plink_files/mass.eigenval", sep=" ", header=F)
pca.mass.diag <- read.table("plink_files/mass.rel.diag", sep=" ", header=F)

for (i in 1:20) {
lmer.all.mass <- lmer.all.mass %>% mutate(!!as.name(paste0("PCA",i,".cor")) :=  !!as.name(paste0("PCA",i)) * (pca.mass.vals[i,]/sum(as.numeric(pca.mass.vals[1:1,1]), na.rm = TRUE)))
lmer.all.mass.variance <- lmer.all.mass.variance %>% mutate(!!as.name(paste0("PCA",i,".var")) :=  !!as.name(paste0("PCA",i)) * (pca.mass.vals[i,]/sum(as.numeric(pca.mass.vals[1:20,1]), na.rm = TRUE)))
}

lmer.all.mass <- lmer.all.mass %>% mutate( PCAallmass.cor = rowSums(.[22:22]))
lmer.all.mass.variance <- lmer.all.mass.variance %>% mutate( PCAallmass.var = rowSums(.[22:22]))
</code></pre>

We can make a summary of how much genetic variance is used over all analysed axis.

<pre class="r"><code>sum(as.numeric(pca.mass.vals[1:1,1]))/sum(as.numeric(pca.mass.diag[,1]))
</code></pre>

Now combine and format the data for the heatmap plot of correlations.

<pre class="r"><code>lmer.all <- data.frame(lmer.all.mass$PCAallmass.cor,lmer.all.tarsus$PCAalltarsus.cor,lmer.all.head$PCAallhead.cor,lmer.all.beak$PCAallbeak.cor,lmer.all.wing$PCAallwing.cor,lmer.all.spleen$PCAallspleen.cor,lmer.all.heart$PCAallheart.cor)
library(data.table)
lmer.all2 <- setDT(lmer.all, keep.rownames = TRUE)[]
lmer.all2[,1] <- c("bio01","bio03","bio04","bio08","bio09","bio12","bio14","bio15","bio18","bio19","elev","mnNDVI")
lmer.all2.long <- melt(setDT(lmer.all2), id.vars = c("rn"), variable.name = "pheno")
</code></pre>

Now plot the heapmap.

<pre class="r"><code>library(ggplot2)
library(viridis)

xlabs.phenos <- c("Mass-SNPs", "Tarsus-SNPs", "Head-SNPs", "Beak-SNPs", "Wing-SNPs", "Spleen-SNPs","Heart-SNPs")

png("Sv6_geno_env_heatpanels_variableaxis.png", width = 800, height = 500)
ggplot(lmer.all2.long, aes(rn, pheno, fill= value)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    scale_y_discrete(labels= xlabs.phenos) +
    theme_minimal(base_size = 18) + 
    theme(axis.text.x = element_text(angle =45, hjust=1),axis.title.x = element_blank(),axis.title.y = element_blank() )
dev.off()
</code></pre>

![screenshot](/Sv6_vignette/Sv6_heatmap.png)


