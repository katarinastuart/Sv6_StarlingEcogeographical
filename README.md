# Sv6_StarlingEcogeographical

This repository contains code related to the following publication:

Stuart KC, Sherwin WB, Cardilini APA, Rollins LA **2022**. Genetics and Plasticity Are Responsible for Ecogeographical Patterns in a Recent Invasion. *Frontiers in Genetics*, 13: 824424, [10.3389/fgene.2022.824424](10.3389/fgene.2022.824424)

In the raw scripts folder you will find:
* 2022_02_12_notebook_16949.pdf: GCTA analysis
* 2022_02_12_notebook_16923.pdf: Morphology-associated alleles correlation analysis (also covered below in the vignette)
* 2022_02_12_notebook_16913.pdf: GradientForest analysis

<h2>Morphology-associated alleles correlation analysis</h2>

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

Here is the code for the anacovis plot.

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





