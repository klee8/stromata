# grab antismash data from genbank file outputs

mkdir results


# telomere to telomere genome

cd Eel728_ONP_polished/
rm ../results/E.elymi_728_ONP_polished_antismash_sum.txt
echo -e "Chr\tregion\treg_start\treg_end\tgene_name\texons\tgene_start\tgene_end\tdesc\ttranslation" > ../results/E.elymi_728_ONP_polished_antismash_sum.txt
for i in *region00*.gbk; do perl ../get_genbank_data.pl $i >> ../results/E.elymi_728_ONP_polished_antismash_sum.txt; done
cd ..


# genomes that the gene_sets were generated from

cd E.festucae_E2368/
rm ../results/E.festucae_E2368_antismash_sum.txt
echo -e "Chr\tregion\treg_start\treg_end\tgene_name\texons\tgene_start\tgene_end\tdesc\ttranslation" > ../results/E.festucae_E2368_antismash_sum.txt
for i in *region00*.gbk; do perl ../get_genbank_data.pl $i >> ../results/E.festucae_E2368_antismash_sum.txt; done
cd ..

cd E.typhina_E8/
rm ../results/E.typhina_E8_antismash_sum.txt
echo -e "Chr\tregion\treg_start\treg_end\tgene_name\texons\tgene_start\tgene_end\tdesc\ttranslation" > ../results/E.typhina_E8_antismash_sum.txt
for i in *region00*.gbk; do perl ../get_genbank_data.pl $i >> ../results/E.typhina_E8_antismash_sum.txt; done
cd ..

cd Eel_masked/
rm ../results/Eel_masked_antismash_sum.txt
echo -e "Chr\tregion\treg_start\treg_end\tgene_name\texons\tgene_start\tgene_end\tdesc\ttranslation" > ../results/Eel_masked_antismash_sum.txt
for i in *region00*.gbk; do perl ../get_genbank_data.pl $i >> ../results/Eel_masked_antismash_sum.txt; done
cd ..

