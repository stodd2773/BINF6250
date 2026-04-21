# Introduction
A profile HMM is a position-specific HMM architecture tailored to model a conserved motif or domain from an MSA. Each alignment column becomes a block of three states (Match *Mi*, Insertion *Ii*, Deletion *Di*), with transitions constrained to flow left-to-right through the motif. This architecture allows the model to capture both column-wise conservation and localized insertions and deletions in homologous sequences.

## Key Features
Position-specific emissions:

- Match states *Mi* use column-specific amino acid distributions estimated from the MSA.

- Insertion states *Ii* share a background composition independent of position.

Explicit gap handling:

- Deletion states *Di* are silent and model gaps by skipping columns without emitting residues.

Left-to-right topology:

Structured transitions between *Mi*, *Ii*, and *Di* encode allowed motif-length variation.

Probabilistic scoring:

- Forward, Viterbi, and Forward–Backward algorithms from `HMM.py` are reused to score sequences and, if desired, retrain parameters.


# Pseudocode

```
profile HMM -> inherit HMM class

    Def __init__(self)

        # alphabet <- letter representations of amino acids
        # pseudocount <- 1
        # self.emit_probs = {}    
        # self.trans_probs = {}     
        # self.begin_probs = {}

    FUNCTION read in fasta for MSA

    FUNCTION assign state to each postion column (match or insertion)

        match if >50% aligned
        insertion if <50% aligned

    FUNCTION label each position in each sequence with state (match, deletion, insertion)

        match or insertion if there is a residue present
            match if global assignment is match
            insertion if global assignment is insertion

        deletion if '-' or symbol denoting gap


    FUNCTION create background distribution (input MSA) # could be in init of profile HMM

        take all observed residues and calc prop for each residue

            make empty background_dist_dict

            for each residue: 

                background_dist = count of (residue + pseudocount) / (total num of residues observed + 20 * pseudocount)
                
                bacground_dist_dict[residue] = background_dist

        return background_dist_dict


    FUNCTION compute emissions (MSA)

    

        for state_index, residue_index in MSA state assignment table.enumerate()
            
        
            if state assignment is a match, execute build_emissions_from_msa for that column using only aligned residues (skipping gaps)
                residues = [seq[residue_index] for seq in msa 
                if seq[residue_index] != '-']

                def build_emissions(self, residues, state_index, emissions):
                        
                    # Count frequencies
                    counts = {}
                    for r in self.alphabet:
                        counts[r] = residues.count(r)
                    total = sum(counts.values())
                    
                    # Compute emission probabilities for match state
                    self.emissions[f'M{state_index}'] = {
                        r: counts[r] + self.pseudocount / total + 20 * self.pseudocount 
                        for r in self.alphabet
                    }
                    
                    # Insertion state uses background distribution
                    emissions[f'I{state_index}'] = self.background_dist()
                
                    return self.emissions


                if state assignment is not match, skip


    FUNCTION calculate transition prob

        initialize with transition_counts with valid transitions: 
            M_i -> M_i+1
            M_i -> I_i 
            M_i -> D_i+1
            ...

        for each labeled state in labeled_state_matrix.shape()
            curr_state <- current labeled state
            next_state <- next state in labeled state matrix

            if curr_state to next_state transition is in transition count dictionary: 
                increment counter for that transition by 1

        
        for each key in the dictionary (original labeled state): 
            sum all counts to get total
            for each key (next state): 
                value (transition count) + psuedocount / total counts + 20 * psuedocount

        
    FUNCTION evaluate (sequences)

        for each seq in sequences: 
            prob = pHMM.forward() to get overall prob
            path = pHMM.viterbi() to get optimal state path for M, D, I

            print prob and path variables 
```

# Successes
The main focus of the group for this project was planning the functions for implementing a profile HMM by extending an existing HMM class. We did not get to fully implementing the code, but in ```HMM.py``` it is clear we have planned each potential function and its output. Our first few meetings entailed discussing and making sure we understood the concept of profile HMMs. We all went through the html provided on the canvas page, walking through each step together, discussing what is being done, and what potential code would look like. We approached the algorithm by each step, planning and writing pseudocode for each. Although we did not finish all of the code, I think what we did accomplish, and the way we went about the project, significantly boosted our understanding of the concept 

# Struggles
This project challenged us conceptually, requiring us to understand and extend code that we did not write. A significant portion of our planning was discussing the concepts and navigating the html link on the canvas page. We initially struggled with understanding what states to assign where, what defines a deletion vs insertion, but open dialogue and using eachother as resources helped us over this hump. 

# Personal Reflections
## Group Leader
Spencer: I was challenged for this project, struggling to understand the difference between assigning residue-specific states and position-specific states, and using this step to inform how to make the profile HMM. I think the planning aspectc of this project went well, and we worked really well together as a team. We made sure everyone was on board conceptually before moving on to pseudocode and writing out the functions for each step. 

## Other member
Sneha: This project was definitely challenging and there were multiple components that were a bit difficult for me to understand. Specifically, it took me a while to understand when insertions and deletions happen and at what position. Another challenge for me was going through the codebase and making sure that I understood everything that was done in order to plan out the psuedocode for the profile HMM accurately. It was nice being able to go through all the resources as a group and talking through each of the steps and having an open dialogue about what we understood and what was confusing.

Marcos: This project was conceptually challenging because of the nature of Profile HMMs. It took me a while to wrap my head around the position-specific states, as well as understanding how to label each position. The potential self-insertion loop, where you could have multiple insertions at the same position was not very intuitive, and took seeing several examples to comprehend. Additionally, we were given an entire codebase to start and work with, which was something different from previous assignments, where we either had no code at all, or had some guiding functions to fill in. In this case, we had to base our subsequent work on the existing code, which meant we first needed to go through and understand what it was doing. We spent all our time on the planning part, wanting to make sure we had a good udnerstanding of pHMMs and how they worked. 

# Generative AI Appendix
As per the syllabus
