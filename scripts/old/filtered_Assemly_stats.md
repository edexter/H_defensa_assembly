### Calculate ASM assembly stats 

````
# Activate environment
conda activate eukaryotic_genome_assembly

# Navigate to project folder
cd aphids

# Calculate ASM stats across all filtered assemblies
while IFS= read -r SAMP; do
  echo "Processing sample ID: $SAMP"

  # Run blobtools filter with the current sample ID
  scripts/asmstats assemblies_filtered/"$SAMP".p_ctg.filtered.fa > assemblies_stats/"$SAMP".stats

done < scripts/sampleID.txt
````

### Aggregate ASM assembly stats

````
RESULTS_DIR="assemblies_stats"
OUTPUT_FILE="assemblies_stats/merged_asm_read_stats.csv"

# Find all files ending in "asm_read_stats" and save them to a temporary file
find "$RESULTS_DIR" -name "*.stats" -print > temp_results.txt

# Add a header to the output file
echo "Filename,sum,n,mean,largest,smallest,N50,L50" > "$OUTPUT_FILE"

# Process each file to extract the required fields and append to the output file
while read -r file; do
  sum=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $2); print $2}' "$file")
  n=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $4); print $4}' "$file")
  mean=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $6); print $6}' "$file")
  largest=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $8); print $8}' "$file")
  smallest=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $10); print $10}' "$file")
  N50=$(awk -F '[=,]' '/N50 =/ {gsub(/ /, "", $2); print $2}' "$file")
  L50=$(awk -F '[=,]' '/N50 =/ {gsub(/ /, "", $4); print $4}' "$file")

  echo "$file,$sum,$n,$mean,$largest,$smallest,$N50,$L50" >> "$OUTPUT_FILE"
done < temp_results.txt

# Clean up temporary file
rm temp_results.txt

echo "File processing complete. Output saved to $OUTPUT_FILE"
````

