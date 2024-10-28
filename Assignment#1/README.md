# bioinformatics-lab
## Bioinformatics Course, First Lab Assignment: Implementation of Sequence Alignment Algorithms
_In this assignment, I will be implementing Needleman-Wunsch and Smith-Waterman algorithms used for Sequence Alignment._

### Sequence Alignment
Sequence alignment is arrangement of DNA, RNA, or protein sequences in order to make it easier to spot the similarities. A similarity between regions hints that there is a functional, structural, or an evolutionary similarity as well. This can mean that those sequences share the same ancestors. On this condition, it can be said that, mismatches indicate point mutations and gaps indicate insertion or deletion mutations. It is possible for us to align sequences which are relatively short, like the one we worked on, during our class. But most of the time, much longer and complicated sequence alignment is needed. To overcome this, algorithms can be used. 

### Needle-man Wunsch Algorithm
This is an global alignment algorithm developed by Saul B. Needleman and Christian D. Wunsch and was published in 1970. During global alignment, every letter in the sequence is aligned, therefore it is needed to look at the entire sequence. 

### Smith-Waterman Algorithm
This is an local alignment algorithm proposed by Temple F. Smith and Michael S. Waterman in 1981. Algorithm compares segments of all possible lengths.