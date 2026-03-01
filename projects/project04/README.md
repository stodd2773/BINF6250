# Introduction

For this project, we implemented a primary assembly algorithm for short-read data: DeBruijn Graphing. This was implemented in simple form where we were to assume perfect sequencing, meaning that everything was sequenced exactly once and there were no errors or variants in the data. 

A DeBruijn graph is composed of nodes and edges, and we used a defaultdict to track edges between the nodes in our graph. This data structure held all edges where the right nodes connected to a left node, stored in a list for that node. We defined the basic class structure for the DeBruijn setup, including descriptions of functions to add and remove edges from the graph.

The purpose of creating our De Bruijn graph was to output a valid sequence from the assembly. This is implemented as a recursive algorithm by considering all valid edges. In a more complex implementation of a Eulerian walk there are heuristics and defined rules for determining the validity of traversing a specific edge in the graph to result in a full graph-traversal. One of these methods is to traverse the graph in a depth first manner to avoid sectioning off any part of the graph in the traversal. 

This project utilizes a FASTQ file with 10 million "perfect" 150 bp reads simulated
from the mouse genome (`GRCm39`).

Assumptions:
Reads interpreted from left to right only (directionality)
No divergent overlaps
No errors
High coverage (coverage = number of reads contributing to overlap)
All reads/sequences are distinct
No gaps or merges at center of read with unmatched bases remaining at end(s)
k is min overlap length required for merge, each read has 2 kmers: left end of read from 0 through k, right end of read from len(read) - k through len(read)
n is actual overlap length between reads, and n >= k
The larger n is to k, the more unique the fit (greater specificity)
At k ~ 15, essentially only one read possible to merge with current (highest specificity)
k will be set as a parameter with default value k=10
Each k-mer has 2 nodes: node, node length = k-1 (left node is 0 through k-1 of kmer, right node is 1 through k of kmer)
4 total nodes per read
Not all reads will merge into a single contig of len(genome)
Attempt to generate longest contig possible by merging reads into single sequence according to overlap length k
Once no further reads can be merged into current contig, new contig is started
It is possible to have unused reads and no possible merges remaining
Once merge is identified, remove merged reads from pool of potential merges
After all reads merged or determined unusable, contigs are compared for potential merges

# Pseudocode
Put pseudocode in this box:

```
1. With gzip open fasta file  in read mode and append the line containing reads to the sequence 
list stripping white spaces
2. Store the k value and initialize self.graph and missing keys get an empty list from default dict
3. Build De Bruijn graph from reads passing in reads and k
4. For every reads skip if shorter than k slide a window across reading one base at a time
5. Split the window into a source node prefix and a destination node suffix and call add edge with source and destination
6. Start eulerian walk by setting a random seed if seed is provided for edge selection 
7. Create an empty path list and push the starting node to the stack
8. Start the iteration from the current node in the stack, if it has outgoing edges, seed is provided and there are multiple outgoing edges to choose from randomly pick the next node
9. If there is only one outgoing edge take the last edge in the list 
10. Remove the edge from the list so it is not visited again and push it onto the stack to continue the iteration for eulerian walk
11. If the current node has no outgoing edges in the stack, dead end is reached and pop off the stack and append the stack to the path 
12. Reverse the path to get the first node order as nodes were added first in last out to the stack and return the path
```

# Successes
As a group we made sure we gave one another time to dicuss the algorithm and psudeocode even if we did not leave with answers at the end. 
We really collerborated well all together in terms of implementation as we improved our code thorughout the meetings we had. Moreover, we helped each other by sharing resources we can get our hands on to understand how the graph alogrithm worked. 


# Struggles
Struggles
We struggled in implementation mainly becuase we were unsure of how edges and nodes are going to be represented. When discussing the psuedocode, this was the biggesest factor of why we could not move forward. We wanted to visualize the alogrithm for easier implementation as well as understanding what we had to work with. This led us to become less aware of our progress and under estimating the remaining time we had left to finish the project.

# Personal Reflections
## Group Leader
Spencer: I struggled the most with this weeks project, with understanding what we wanted our code to do. I personally struggled with converting a visual interpretation of a graph to an executable program that assembles reads. I think that a huge success in the face of these challenges was our groups ability to collaborate, share resources, talk through concepts and lift each other up. Through many meetings and collaboration we were able to work out a solution and implement 

## Other members
Thu Thu Han - I personally believe that this week's project was one of the most challenging ones. I struggled during the implementation of the nodes and edges and part of it is becuase I struggled with visualizing the graph. As a group we went through many resources to understand the DeBruijn Graph Algorithm and I think the resources my teammate has shared really helped me to start understanding. I think our meeting as a group with Marcus has helped a lot as well as he gave us hints of how to implement the graph and gave us a brief overview of the struggles we should expect from implementation. I think with this project I really learned alot in many aspects. I got to explore new data structures like stack and new algorithms like the sliding window algorithm.

Stefanie Moreno- This week's project pushed us to integrate several core concepts in genome assembly while also managing the practical challenges of working with a high-volume sequencing dataset. One of the biggest hurdles was designing a pipeline that could handle millions of reads without exhausting memory, which required transitioning from a batch‑loading approach to a streaming FASTQ reader and carefully building the De Bruijn graph one read at a time. Debugging the Eulerian traversal and ensuring that connected components were correctly identified also took time, especially because small logic errors could silently break the assembly. Despite these challenges, it was rewarding to see the full pipeline come together: the graph construction, Eulerian walks, contig generation, and output files all worked cohesively. The project helped solidify my understanding of how k‑mers, graph connectivity, and Eulerian paths interact to reconstruct genomic sequence, and it gave me a much deeper appreciation for the complexity of real genome assembly tools. 

Our team met many times and we worked really well together- Spencer has been group leader for each project this semester, so he was very adept at setting up our repo, organizing our meetings, and ensuring great collaboration skills from all of us. I worked with Thu Thu on Project02 and will be working with her again for Project05- she is a highly intelligent student and is great at explaining things when I have a question or don't quite understand something. We struggled at first with identifying what was to be nodes and edges, but after meeting as a team with Marcus during office hours, we felt much more confident in our understanding of how to approach the assembly. We intended on implementing a few other steps in our script, but we were short on time and had to settle with what we had. 

We decided to output three files- one FASTA file containing the contigs, a csv file with the metadata for each contig: contig ID, start position, and coverage depth, and finally, and an easy to read summary stats text file that contained input file name, kmer size, random seed, number of reads processed, average read length, estimated coverage graph statistics (nodes, edges), and assembly results (number of contigs, total assembly length, N50, assembly fraction and top 20 contig lengths). Our reason for creating these output files was because most of the results were to be printed to the terminal per the driver program, and since our peer reviewers are not running our code, we decided to share the details of our assembly in file format.


# Generative AI Appendix
As per the syllabus
