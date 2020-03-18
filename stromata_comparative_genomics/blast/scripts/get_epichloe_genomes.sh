# retrieve fasta files for Epichloe genomes from NCBI (note search for the species name via NCBI genome database)

wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/223/075/GCA_000223075.2_EpiAma_1.0/GCA_000223075.2_EpiAma_1.0_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/877/375/GCA_000877375.1_EamaE4668v1/GCA_000877375.1_EamaE4668v1_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/729/855/GCA_000729855.1_ASM72985v1/GCA_000729855.1_ASM72985v1_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/729/845/GCA_000729845.1_EbacATCC_200745v1/GCA_000729845.1_EbacATCC_200745v1_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/222/915/GCA_000222915.1_EpiBra_1.0/GCA_000222915.1_EpiBra_1.0_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/008/055/GCA_001008055.1_ASM100805v1/GCA_001008055.1_ASM100805v1_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/008/065/GCA_001008065.1_ASM100806v1/GCA_001008065.1_ASM100806v1_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/729/905/GCA_000729905.1_EbroATCC_200750v1/GCA_000729905.1_EbroATCC_200750v1_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/315/335/GCA_000315335.1_EElymi1.0/GCA_000315335.1_EElymi1.0_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/591/845/GCA_002591845.1_ASM259184v1/GCA_002591845.1_ASM259184v1_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/151/575/GCA_000151575.2_ASM15157v2/GCA_000151575.2_ASM15157v2_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/814/445/GCA_003814445.1_ASM381444v1/GCA_003814445.1_ASM381444v1_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/163/435/GCA_002163435.1_ASM216343v1/GCA_002163435.1_ASM216343v1_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/222/895/GCA_000222895.2_NeGans_1.0/GCA_000222895.2_NeGans_1.0_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AF/RG/AFRG01/AFRG01.1.fsa_nt.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AM/DK/AMDK01/AMDK01.1.fsa_nt.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/JF/GW/JFGW01/JFGW01.1.fsa_nt.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LC/TT/LCTT01/LCTT01.1.fsa_nt.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/308/955/GCA_000308955.1_ETyphina_1.0/GCA_000308955.1_ETyphina_1.0_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/222/955/GCA_000222955.2_AciTa_1.0/GCA_000222955.2_AciTa_1.0_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LE/LE/LELE01/LELE01.1.fsa_nt.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/318/955/GCA_002318955.1_ASM231895v1/GCA_002318955.1_ASM231895v1_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/226/195/GCA_000226195.2_EfestFl1_2.0/GCA_000226195.2_EfestFl1_2.0_genomic.fna.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/319/005/GCA_002319005.1_ASM231900v1/GCA_002319005.1_ASM231900v1_genomic.fna.gz"


mkdir data
mv	GCA_000223075.2_EpiAma_1.0_genomic.fna	data/E.amarillans_E57.fa
mv	GCA_000877375.1_EamaE4668v1_genomic.fna	data/E.amarillans_E4668.fa
mv	GCA_000729855.1_ASM72985v1_genomic.fna	data/E.aotearoae_E899.fa
mv	GCA_000729845.1_EbacATCC_200745v1_genomic.fna	data/E.baconii_E1031.fa
mv	GCA_000222915.1_EpiBra_1.0_genomic.fna	data/E.brachyelytri_E4804
mv	GCA_001008055.1_ASM100805v1_genomic.fna	data/E.bromicola_E7561.fa
mv	GCA_001008065.1_ASM100806v1_genomic.fna	data/E.bromicola_AL0434.fa
mv	GCA_000729905.1_EbroATCC_200750v1_genomic.fna	data/E.bromicola_E502.fa
mv	GCA_000315335.1_EElymi1.0_genomic.fna	data/E.elymi_E56.fa
mv	GCA_002591845.1_ASM259184v1_genomic.fna	data/E.elymi_E757.fa
mv	GCA_000151575.2_ASM15157v2_genomic.fna	data/E.festucae_E2368.fa
mv	GCA_003814445.1_ASM381444v1_genomic.fna	data/E.festucae_E894.fa
mv	GCA_002163435.1_ASM216343v1_genomic.fna	data/E.festucae_Fg1.fa
mv	GCA_000222895.2_NeGans_1.0_genomic.fna	data/E.gansuensis_E7080.fa
mv	AFRG01.1.fsa_nt	data/E.glyceriae_E277.fa
mv	AMDK01.1.fsa_nt	data/E.inebrians_E818.fa
gunzip  JFGW01.1.fsa_nt.gz
mv	JFGW01.1.fsa_nt	data/E.mollis_E3601.fa
mv	LCTT01.1.fsa_nt	data/E.sylvatica_E7368.fa
mv	GCA_000308955.1_ETyphina_1.0_genomic.fna	data/E.typhina_E8.fa
mv	GCA_000222955.2_AciTa_1.0_genomic.fna	data/E.typhina_spp.Poae_E5819.fa
mv	LELE01.1.fsa_nt	data/E.uncinata_e167.fa
mv	GCA_002318955.1_ASM231895v1_genomic.fna	data/E.festucae_AL9436.fa
mv	GCA_000226195.2_EfestFl1_2.0_genomic.fna	data/E.festucae_Fl1_2019.fa
mv	GCA_002319005.1_ASM231900v1_genomic.fna	data/E.bromicola_E7626.fa


##### NOTE : did a sanity check of number of contigs for all files, two did not match NCBI metadata
# E. mollis	E3601	AL9924	sexual	Schardl et al. 2014	Holcus mollis	PRJNA215230	SAMN02911893		(27871 contigs instead of metadata 1068)
# E. festucae	Fl1 (2019)					PRJNA376182	SAMN02981348 	35.1499	 (759 contigs instead of metadata 1398)



##### NOTE : background on genomes
#E. amarillans		E57	sexual		Schardl et al. 2013		PRJNA67301 AFRF00000000.1
#wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AF/RF/AFRF01/AFRF01.1.fsa_nt.gz"
#		E4668	sexual		Kentucky data base JFGZ00000000.1
#wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/JF/GZ/JFGZ01/JFGZ01.1.fsa_nt.gz"

#E. aotearoae		E899	asexual		Kentucky data base

#E. baconii		E1031	sexual		Kentucky data base JFGY00000000.1
#wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/JF/GY/JFGY01/JFGY01.1.fsa_nt.gz"

#E. brachyelytri		E4804	sexual		Schardl et al. 2013		PRJNA67245
#wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AF/RB/AFRB01/AFRB01.1.fsa_nt.gz"

#E. bromicola		E7561/AL0426/2	sexual		Kentucky database PRJNA274992
#wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LB/NH/LBNH01/LBNH01.1.fsa_nt.gz"
#		ALO434	sexual		Kentucky database  	PRJNA274994
#wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LB/NI/LBNI01/LBNI01.1.fsa_nt.gz"
#		E502	sexual		Kentucky database PRJNA221343 JFHA00000000.1 Epichloe bromicola ATCC 200750
#wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/JF/HA/JFHA01/JFHA01.1.fsa_nt.gz"

#E. elymi		E56	sexual		Schardl et al. 2013		PRJNA173776
#wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AM/DJ/AMDJ01/AMDJ01.1.fsa_nt.gz"
#		E757	sexual		Scott unpublished		SAMN07182378                        
#		E6255/NFE728	sexual                                                          <<<<<<<<<<< Don't have this one

#E. festucae		E2368	sexual		Schardl et al. 2013		PRJNA42133
#wget "ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AD/FL/ADFL02/ADFL02.1.fsa_nt.gz"
#		Fl1	sexual		Schardl et al. 2013		PRJNA51625
#		Fg1	sexual		Scott unpublished		SAMN06710538

#E. festucae var lolii		AR37	asexual		Agresearch                              <<<<<<<<<<< Get from Agresearch                      
#		AR5	asexual		Agresearch
#		AR48	asexual		AgResearch
#		AR1	asexual		AgResearch

#E. glyceriae		E277	sexual		Schardl et al. 2013		PRJNA67247
#		E2772	sexual		Schardl et al. 2013

#E. mollis		E3601	sexual		Kentucky database

#E. sylvatica		E7368	sexual		Kentucky database

#E. typhina		E8	sexual		Schardl et al. 2013		PRJNA174036

#E. typhina spp. Poae		E5819	sexual		Schardl et al. 2013		PRJNA68441

#E. inebrians		E818	asexual		Schardl et al. 2013		PRJNA174039
#E. gansuensis		E7080	asexual		Schardl et al. 2013		PRJNA67299
