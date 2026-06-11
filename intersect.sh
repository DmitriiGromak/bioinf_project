GFF="genomic.gff"
GENOME_SIZES="genome.sizes"
PROM=1000
DOWN=200
awk '/^[^#]/{print $1"\t"$5}' "$GFF" \
    | grep "^NC_" \
    | sort -k1,1 -u > "$GENOME_SIZES"
awk '$3=="exon" && $1~/^NC_/' "$GFF" \
    | awk '{OFS="\t"; print $1,$4-1,$5}' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > exons.bed
awk '$3=="gene" && $1~/^NC_/' "$GFF" \
    | awk '{OFS="\t"; print $1,$4-1,$5}' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    | bedtools subtract -a - -b exons.bed > introns.bed
awk '$3=="gene" && $1~/^NC_/' "$GFF" \
    | awk 'NR==FNR{sizes[$1]=$2; next}
           {OFS="\t"; strand=$7; s=$4-1; e=$5; chr=$1; L=sizes[chr]
            if(strand=="+"){ps=(s-1000>0?s-1000:0); pe=s}
            else            {ps=e; pe=(e+1000<L?e+1000:L)}
            if(pe>ps) print chr,ps,pe}' \
      "$GENOME_SIZES" - \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > promoters.bed
awk '$3=="gene" && $1~/^NC_/' "$GFF" \
    | awk 'NR==FNR{sizes[$1]=$2; next}
           {OFS="\t"; strand=$7; s=$4-1; e=$5; chr=$1; L=sizes[chr]
            if(strand=="+"){ps=e; pe=(e+200<L?e+200:L)}
            else            {ps=(s-200>0?s-200:0); pe=s}
            if(pe>ps) print chr,ps,pe}' \
      "$GENOME_SIZES" - \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > downstream.bed
cat exons.bed introns.bed promoters.bed downstream.bed \
    | cut -f1-3 \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    | bedtools complement -i - -g "$GENOME_SIZES" > intergenic.bed
for src in zhunt_filtered.bed zdna_bomhunt.bed g4_predictions.bed; do
    tag="${src%%.*}"
    grep "^NC_" "$src" | cut -f1-3 | awk '$3>$2' \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - > "${tag}_merged.bed"
    echo "$tag: $(wc -l < $src) raw → $(wc -l < ${tag}_merged.bed) merged"
done
cp zhunt_filtered_merged.bed zhunt_merged.bed
cp zdna_bomhunt_merged.bed   zdnabert_merged.bed
cp g4_predictions_merged.bed g4_merged.bed
echo -e "Feature\tZhunt_count\tZDNABERT_count\tG4_count" > table1_counts.txt
for FEAT in exons introns promoters downstream intergenic; do
    Z=$(bedtools intersect -a zhunt_merged.bed    -b ${FEAT}.bed -u | wc -l)
    D=$(bedtools intersect -a zdnabert_merged.bed -b ${FEAT}.bed -u | wc -l)
    G=$(bedtools intersect -a g4_merged.bed       -b ${FEAT}.bed -u | wc -l)
    echo -e "${FEAT}\t${Z}\t${D}\t${G}" >> table1_counts.txt
done
echo -e "Feature\tTotal_intervals\tZhunt_with_hit\tZDNABERT_with_hit\tG4_with_hit" > table2_counts.txt
for FEAT in exons introns promoters downstream intergenic; do
    TOTAL=$(wc -l < ${FEAT}.bed)
    Z=$(bedtools intersect -a ${FEAT}.bed -b zhunt_merged.bed    -u | wc -l)
    D=$(bedtools intersect -a ${FEAT}.bed -b zdnabert_merged.bed -u | wc -l)
    G=$(bedtools intersect -a ${FEAT}.bed -b g4_merged.bed       -u | wc -l)
    echo -e "${FEAT}\t${TOTAL}\t${Z}\t${D}\t${G}" >> table2_counts.txt
done