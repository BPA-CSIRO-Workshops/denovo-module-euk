#Commands for DeNovo Training
***

#Exercise #1: Fusarium first pass with a goal

##Goal: Identify a fusarium sample is "closer" to F. graminearum or F.pseudograminearum

Previous knowledge:

* F. graminareum has a cluster producing PKS6 and NRSP7, while F. pseudograminareum produces PKS40 and NRPS32

Data:

* Proteins sequences for:
 * F. graminareum (non necrotrophic): PKS6 and NRSP7
  * F. pseudograminareum (necrotrophic): PKS40 and NRPS32
* Blast database of cereal pathogen proteins.

Strategy:

* Check PKS6-NRSP7 and PKS40-NRPS32 cluster presence.

Assembly goals:

* Assembly goal (I): to capture a good enough representation of the protein-coding space to get blast matches
* Assembly goal (II): to accurately represent the relevant whole cluster loci in a single sequence.

##Task 1.1: First pass assembly, k=71

First, assembly:

~~~
cd denovo/fusarium
mkdir abyss_k71
cd abyss_k71
abyss-pe in="../CS3270_A8733_GCCAAT_L001_R1.fastq ../CS3270_A8733_GCCAAT_L001_R2.fastq" k=71 name=CS3270_abyss_k71 np=4 > CS3270_abyss_k71.log 2>&1
~~~

Stats... Ok, no stats, but we can always use abyss-fac

~~~
abyss-fac CS3270_abyss_k71-contigs.fa
~~~
n	|n:500	|L50	|min	|N80	|N50	|N20	|E-size	|max	|sum	|name
------	|------	|------	|------	|------	|------	|------	|------	|------	|------	|------
27	|13	|2	|970	|6004	|13202	|52602	|28712	|52602	|112849	|CS3270\_abyss\_k71-unitigs.fa
5	|1	|1	|128429	|128429	|128429	|128429	|128429	|128429	|128429	|CS3270\_abyss\_k71-contigs.fa

Strange! Time for some analysis:

- Check frequencies for kmers kept/discarded/etc.
- Check spectra-cn and compare with expectations.

~~~
less CS3270_abyss_k71.log 
less coverage.hist

gnuplot
gnuplot> set xrange [0:50]
gnuplot> set yrange [0:4000000]
gnuplot> plot "coverage.hist"

gnuplot
gnuplot> set xrange [0:200]
gnuplot> set yrange [0:5000]
gnuplot> plot "coverage.hist"

kat comp -o reads_vs_abyss_k71 -t 4 -C -D '../*.fastq' CS3270_abyss_k71-contigs.fa
~~~
Looks like we are not assembling this bit, let's have another look at the spectra
~~~
kat comp -o reads_vs_abyss_k71 -t 4 -C -D --d1_bins 2000 '../*.fastq' CS3270_abyss_k71-contigs.fa
kat plot spectra-cn -y 600 -x 2000 -o reads_vs_abyss1-main.mx.spectra-cn_2000.png reads_vs_abyss_k71-main.mx
~~~

Take the output and blast it in NCBI. What is it? Surprising?

This assembly will get us nowhere, let's choose a lower K to gain coverage and start again.

##Task 1.2: First pass assembly, k=27
~~~
cd denovo/fusarium
mkdir abyss_k27
cd abyss_k27
abyss-pe in="../CS3270_A8733_GCCAAT_L001_R1.fastq ../CS3270_A8733_GCCAAT_L001_R2.fastq" k=27 name=CS3270_abyss_k27 np=4 > CS3270_abyss_k27.log 2>&1
~~~

Stats look better:

n      |n:500  |L50  |min  |N80     |N50     |N20      |E-size  |max      |sum      |name
---    |---    |---  |---  |---     |---     |---      |---     |---      |---      |---
30645  |2717   |430  |502  |11354   |25336   |47966    |31027   |147694   |36.14e6  |CS3270\_abyss\_k27-unitigs.fa
21511  |350    |33   |527  |157565  |338989  |630228   |407098  |1265237  |36.52e6  |CS3270\_abyss\_k27-contigs.fa
21327  |205    |17   |527  |332444  |716132  |1265237  |791882  |1880850  |36.51e6  |CS3270\_abyss\_k27-scaffolds.fa

Let's check a bit anyway:

~~~
less CS3270_abyss_k27.log
less coverage.hist
gnuplot
gnuplot> set xrange [0:50]
gnuplot> set yrange [0:4000000]
gnuplot> plot "coverage.hist"

kat comp -o reads_vs_abyss_k27 -t 4 -C -D '../*.fastq' CS3270_abyss_k27-scaffolds.fa
~~~

###Question: any tools you can use to check kmer spectra at any K before assembling?

###Question: can you predict what will happen if you use KAT with larger K values?

##Task 1.3: Will the assembly answer the biological question?

Use blast and the proteins to check...

***

#Exercise #2: Chalara scaffolding using LMP
##Task 2.1: let's have a look at the PE assembly
~~~
abyss-pe name=cha1 k=27 in="../LIB2570_raw_R1.fastq ../LIB2570_raw_R2.fastq" np=4
~~~

n       |n:500  |L50   |min  |N80   |N50   |N20    |E-size  |max    |sum      |name
---     |---    |---   |---  |---   |---   |---    |---     |---    |---      |---
394789  |11681  |2094  |500  |2880  |6254  |12303  |7936    |44414  |44.07e6  |cha1-unitigs.fa
394199  |11673  |2097  |500  |2887  |6255  |12303  |7937    |44414  |44.11e6  |cha1-contigs.fa
394161  |11647  |2095  |500  |2898  |6269  |12303  |7944    |44414  |44.12e6  |cha1-scaffolds.fa
~~~
abyss-pe name=cha2 k=61 in="../LIB2570_raw_R1.fastq ../LIB2570_raw_R2.fastq" np=4
~~~
n       |n:500  |L50   |min  |N80   |N50   |N20    |E-size  |max    |sum      |name
---     |---    |---   |---  |---   |---   |---    |---     |---    |---      |---
130547  |12596  |1770  |500  |3352  |8379  |16380  |10518   |54300  |50.75e6  |cha2-unitigs.fa
130210  |12575  |1771  |500  |3363  |8382  |16380  |10525   |54300  |50.78e6  |cha2-contigs.fa
130182  |12555  |1769  |500  |3377  |8394  |16380  |10534   |54300  |50.78e6  |cha2-scaffolds.fa

Contiguity is worse than fusarium, why?

Answers:

- Genome characteristics.
- Data:
 - Coverage
 - Error distributions
 - Read sizes
 - Fragment sizes
 - Kmer spectra


***
##Task 2.2: Let's put some LMP in there.
~~~
abyss-pe name=chalmp1 k=61 in="../LIB2570_raw_R1.fastq ../LIB2570_raw_R2.fastq" mp="lmp1" lmp1="../LIB8209_raw_R1.fastq ../LIB8209_raw_R2.fastq" np=4
~~~
n       |n:500  |L50   |min  |N80   |N50   |N20    |E-size  |max    |sum      |name
---     |---    |---   |---  |---   |---   |---    |---     |---    |---      |---
130548  |12596  |1770  |500  |3352  |8379  |16380  |10518   |54300  |50.75e6  |chalmp1-unitigs.fa
130211  |12575  |1771  |500  |3363  |8382  |16380  |10525   |54300  |50.78e6  |chalmp1-contigs.fa
130148  |12545  |1769  |500  |3380  |8400  |16380  |10535   |54300  |50.78e6  |chalmp1-scaffolds.fa

So... what happened?

Data:

- Kmer spectra
- Fragment sizes
- Any hints on the protocol?

A not-so-obvious property:

~~~
kat plot spectra-mx --intersection -x 30 -y 15000000 -o pe_vs_lmp-main.mx.spectra-mx.png pe_vs_lmp-main.mx
~~~


LMP require pre-processing, right?

***
##Task 2.3: Let's try with processed LMP

Prior task (already made) preprocess the LMP with nextclip.

~~~
abyss-pe name=chalmpproc1 k=61 in="../../cha_raw/LIB2570_raw_R1.fastq ../../cha_raw/LIB2570_raw_R2.fastq" mp="proclmp1" proclmp1="../LIB8209_preproc_R1.fastq ../LIB8209_preproc_R2.fastq" np=4
~~~

n       |n:500  |L50   |min  |N80    |N50    |N20     |E-size  |max     |sum      |name
---     |---    |---   |---  |---    |---    |---     |---     |---     |---      |---
130548  |12596  |1770  |500  |3352   |8379   |16380   |10518   |54300   |50.75e6  |chalmpproc1-unitigs.fa
130211  |12575  |1771  |500  |3363   |8382   |16380   |10525   |54300   |50.78e6  |chalmpproc1-contigs.fa
122061  |6679   |167   |500  |18609  |87510  |187171  |106178  |397967  |51.15e6  |chalmpproc1-scaffolds.fa

That's much better!

***
#Excercise #3: Chalara: beyond first pass

Do you remember these?

- Genome characteristics.
- Data:
 - Coverage
 - Error distributions
 - Read sizes
 - Fragment sizes
 - Kmer spectra

Let's think a bit and try to improve the assembly...

###Question: if you look at the pre-processed LMP, do you notice anything peculiar?


##Example: Testing the inclussion of heavily pre-procesed LMP coverage into the DBG
~~~
abyss-pe name=chalmp2 k=61 se="../LIB8209_preproc_single.fastq" lib="lmp1" lmp1="../LIB8209_preproc_R1.fastq ../LIB8209_preproc_R2.fastq" np=4 >chaproclmp2.log 2>&1
~~~

n       |n:500  |L50   |min  |N80    |N50     |N20     |E-size  |max     |sum      |name
---     |---    |---   |---  |---    |---     |---     |---     |---     |---      |---
128306  |10150  |1182  |500  |4954   |12585   |24777   |15693   |68684   |51.03e6  |chaproclmp4-unitigs.fa
118939  |5772   |362   |500  |15717  |41870   |82033   |52819   |252066  |52.24e6  |chaproclmp4-contigs.fa
116139  |4014   |141   |500  |41145  |109273  |224986  |133450  |460979  |52.56e6  |chaproclmp4-scaffolds.fa

