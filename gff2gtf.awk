#!/usr/local/env awk -f
{
    if(match($0, "ID \\\"([^\\\"]+)\\\";?", id)) {
        if($3 == "gene") {
            last_gene = id[1];
            last_transcript = ""
        };
        if($3 == "transcript" || $3 == "mRNA") {
            last_transcript = id[1]
        };
        for(i=1;i<=NF;i++){
            if(i==9) {
                printf("gene_id \"%s\"; transcript_id \"%s\"; ", last_gene, last_transcript);
            }
            printf($i);
            if(i<9) {
                printf("\t");
            }
            if(i>=9) {
                printf(" ")
            }
        }
        print ""
    }
}
