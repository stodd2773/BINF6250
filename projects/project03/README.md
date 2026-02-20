# Introduction
This project implements a Gibbs Sampling algorithm to identify DNA sequence motifs. The script uses a probabilistic Markov Chain Monte Carlo (MCMC) approach to iteratively find candidate motifs and converge on an enriched sequence pattern represented as a Position Frequency Matrix (PFM).
# Pseudocode
```
Import  libraries

Load BAM file
For each read in BAM file:
    Extract the read sequence
    Store the sequence in list "seqs"

Load genome FASTA file
Load GFF annotation file

Initialize an empty list "seqss"

For each sequence in FASTA:
    For each annotation entry in GFF:
        If annotation type is "CDS":
            Extract promoter sequence (50 bp upstream)
            If promoter contains "AGGAGG":
                Add promoter sequence to "seqss"

Print number of sequences in seqs and seqss

initialize_random_motif(seqs, k):

    Create an empty list motif_sequence_list

    For each sequence in seqs:
        Randomly choose start index between 0 and (length of sequence − k)
        Extract substring of length k starting at random index
        Add substring to motif_sequence_list

    Return motif_sequence_list

Function initialize_random_motif(seqs, k):

    Create empty list motif_sequence_list

    For each sequence in seqs:
        Randomly choose start index between 0 and (length of sequence − k)
        Extract substring of length k starting at random index
        Add substring to motif_sequence_list

    Return motif_sequence_list

select_new_motif(motif_list, seqs, k):

    Randomly choose one sequence index to omit
    Remove its motif from motif_list temporarily

    Build new PFM from remaining motifs
    Build new PWM from new PFM

    Initialize empty list scores

    For each possible k-mer position in omitted sequence:
        Extract forward k-mer
        Compute reverse complement k-mer

        Calculate PWM score for forward k-mer
        Calculate PWM score for reverse complement

        Add both scores to scores list

    Apply softmaxxing:
        For each score:
            Compute 2^score
        Divide each by total sum

    Randomly select new motif index using normalized probabilities

    If selected index corresponds to forward strand:
        Extract k-mer from original sequence
    Else:
        Extract reverse complement k-mer
        Update sequence to reverse complement

    Return omitted sequence index and new motif

GibbsMotifFinder(seqs, k, seed, max_iterations, ic_threshold, convergence_iterations):

    Set random seed

    Initialize motif_list using initialize_random_motif

    Build initial PFM from motif_list
    Build initial PWM
    Compute initial information content (IC)

    Set iteration counter i = 1
    Set convergence_counter = 0

    While i <= max_iterations:

        Call select_new_motif
        Update motif_list with returned motif

        Build new PFM
        Compute new IC

        Calculate absolute difference between new IC and current IC

        If difference < ic_threshold:
            Increment convergence_counter
        Else:
            Reset convergence_counter to 0

        If convergence_counter == convergence_iterations:
            Return new PFM

        Update current IC
        Increment iteration counter

    If max_iterations reached:
        Print "failed to converge"
        Return latest PFM
```

# Successes
Solid teamwork drove some major successes and key takeaways while working on this project. 
Perhaps the biggest driver of success, was the teams ability to come together and effectively plan. In our first meeting, we had a conceptually driven discussion regarding the application of the Gibbs sampling algorithm, familiarizing ourselves with what we wanted our output to look like, and how the Gibbs sampling steps outlined in lecture would take us there. Each step along the way was transcribed into our own words, which created a 'pseudo-psuedocode' outline, from which translating into pseudocode was made simpiler. Taking these fragmented steps in our planning set us up for easier implementation, not having to concern ourselves with conceptual understanding as the code was implemented. 


# Struggles
In terms of the algorithm, the part we struggled with the most was the socring function and dealing with negative numbers for our probabilities. After some discussion and research, we came across the softmax equation, which involved using the power `e` to make sure all the numbers were positive. After this was discussed in class, and we realized we were using `log2` in the scoring function, we changed this to use the power of 2 instead. We also struggled with the data preparation that we worked on during the second week of the assignment. Using `MACS2` peak calling and then querying a reference genome based on the BED file produced was something new to us. We used the tutorial provided during class to familiarize ourselves with the tool, and adapted to our use case. We were able to extract ~500 sequences from the peaks, which we then fed into our algorithm. However, we noticed that the motif the algorithm was identifying was full of `A` and `T` bases, which seemed odd to us. We thought it might have been identyfing poly-A or -T tails, instead of the actual motif. After some tweaking and playing around, we were not able to make it work, so we decided to test our algorithm with the preprocessed CDS data extracted from the .gff file. While we weren't successful, it was a learning experience to familiarize ourselves with a new tool and bioinformatics process.

# Personal Reflections
## Group Leader
Spencer - This project felt conceptually more difficult to grasp than it was to implement. I felt that Marcos, Sneha, and I collaborated particularly well, with key emphasis on our planning meetings which ensured each of us had an understanding of the Gibbs sampling algorithm and how it applied to motif finding. These collaborative meetings, which we were able to have a handful of times, highlighted the efficacy of group work to deepen understanding of complex topics. Each member of the group brought their own understanding and perspective of the project, which was crucial to filling in the gaps of my understanding. This brought about effective planning, which I thought gave way to implementation that was much less difficult than originally anticipated.  

## Other members
Sneha- Spencer and Marcos were both great to work with. We were able to meet a handful of times throughout the two weeks to plan and implement our code rather than using a divide and conquer strategy. The planning phase was especially helpful to make sure we were all on the same page about what we had learned in class and what needed to be done for our project and specific functions.  

Marcos - this was a challenging but very rewarding project. and both Sneha and Spencer were very helpful to make this work as well as it did. As mentioned, we met several times throughout the two weeks, during which we thought out and planned the outline of the algorithm together, as well as doing the coding at the same time. This worked really well, as it allowed all of us to soak in the content and help each other when there were doubts. Once we got that working, we tried to use `MACS2` for peak calling to prepare a more curated dataset. We put in a lot of time to figure out how the tool worked and how to use it for our project, but ended up not having a successful dataset. Despite that, I think it was a very valuable experience to learn how to use a new tool. 

# Generative AI Appendix
As per the syllabus
