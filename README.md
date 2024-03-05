This program was created for an assignment in Bi625 for the Knight Campus Bioinformatics and Genomics Master Program.

It will take in a fasta file and a text file listing one motif per line and return a png graph with the motifs marked on the given sequences.
Currently the code only supports up to 10 motifs at a time. It can take in any number of sequences, but formatting may start having errors for sequence counts greater than 10.

Arguments:
-f, --fasta_file: Path to sorted fasta file
-m, --motifs_file: help="Path to motif file", type=str, required=True)
-t, --theme: Color theme for image. 
            Options are classic, light, dark, mono, contrast, viridi(color blind friendly), and vintage(default).
-s, --stagger: Set to True to stagger motifs below sequence rather than have them overlap the sequence.
-d, --dense: Set to True to have motifs be solid and more transparent, rather than hollow rectangles.
