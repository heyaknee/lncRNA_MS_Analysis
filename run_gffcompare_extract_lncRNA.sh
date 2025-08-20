#set paths
BASE_DIR="path_to_your_base_directory"
ASSEMBLY_DIR="$BASE_DIR/assembly/merged_gtfs"
GENOME_GTF="$BASE_DIR/genome/gencode.v47.annotation.gtf"
OUTPUT_DIR="$BASE_DIR/assembly/gffcompare_output"
NOVEL_DIR="$BASE_DIR/assembly/lncRNA_transcripts"

mkdir -p $OUTPUT_DIR
mkdir -p $NOVEL_DIR

# List all condition folders
CONDITIONS=("Non-MS_control" "Non-MS" "PPMS" "RRMS" "SPMS")

# Loop through each condition and run gffcompare + extract novel lncRNAs
for CONDITION in "${CONDITIONS[@]}"; do
	    COND_DIR="$ASSEMBLY_DIR/$CONDITION"
	        
	        if [ -d "$COND_DIR" ]; then
			        for GTF_FILE in "$COND_DIR"/*.gtf; do
					            SAMPLE_NAME=$(basename "$GTF_FILE" .gtf)
						                echo "🔹 Running gffcompare for $SAMPLE_NAME in $CONDITION..."

								            # Run gffcompare
									                gffcompare -r "$GENOME_GTF" -o "$OUTPUT_DIR/${SAMPLE_NAME}_gffcompare" "$GTF_FILE"

											            # Extract novel lncRNAs (class_code "u")
												                echo "🧬 Extracting novel lncRNAs from $SAMPLE_NAME..."
														            awk '$3 == "transcript" && /class_code "u"/' "$OUTPUT_DIR/${SAMPLE_NAME}_gffcompare.annotated.gtf" > "$NOVEL_DIR/${SAMPLE_NAME}_novel_lncRNAs.gtf"
															            done
																        else
																		        echo "⚠️ Warning: Directory $COND_DIR not found!"
																			    fi
																		    done

																		    echo "✅ GFFCompare & Novel lncRNA Extraction Completed!"

