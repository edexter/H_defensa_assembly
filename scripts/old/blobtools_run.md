











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
  C=$(awk -F '[=,]' '/C =/ {gsub(/ /, "", $2); print $2}' "$file")
  G=$(awk -F '[=,]' '/G =/ {gsub(/ /, "", $2); print $2}' "$file")

  echo "$file,$sum,$n,$mean,$largest,$smallest,$N50,$L50,$C,$G" >> "$OUTPUT_FILE"
done < temp_results.txt

# Clean up temporary file
rm temp_results.txt

echo "File processing complete. Output saved to $OUTPUT_FILE"
````



````
GC: 0.35 - 0.50
Superkingdom != Eukaryota
Family != Erwiniaceae
Family != Yersiniaceae
Family != Morganellaceae
Genus != Providencia
Phylum == Uroviricota
````











