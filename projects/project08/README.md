# Introduction
This project implements the viterbi algorithm given an HMM model. The viterbi algorithm finds the optimal path of hidden states given a sequence of observations. It also implements the forwards, backwards, and forwards-backwards algorithms, which can be used for calculating sequence likelihoods and posterior probabilities. 

# Pseudocode
Pseudocode for Viterbi algorithm:

```

Class HMM(emission_prob, transition_prob, initial_prob):

    FUNCTION viterbi(self, observations):
    
    INITIALIZATION: start with initialization of probabilities for each state for each observation
    initialize two matrices, probability matrix and move matrix 
    first column of probability matrix:
        for each state: 
            initial probability of state k * emission probability for observation 0 in state k

    RECURSION: 
        For each remaining observation: 
            for each possible state:
                previous probability for preceding state  * transition probability * emission probability
                take max value after all states are iterated through
                store state that contributed to max value
                return score matrix, move matrix, max of last column (or Termination function)
                
    TERMINATION:
        after all observations are iterated through
        Take max from last column and begin traceback
    
    TRACEBACK:
        take state contributing to max of last column
        add to state path list 
        move one back in the matrix and look at cell in corresponding state row
        get state contributing to max in new cell
        add to state path list
        move one back in matrix and look at cell in corresponding state row
        repeat for all observations
        reverse list of state path
        return state path list (final output)

    return output from traceback
```

Pseudocode for Forward, Backward, and Forward-Backward algorithms:

```
Class HMM(emission_prob, transition_prob, initial_prob):

    FUNCTION forward(self, observations):

        INITIALIZATION: start with initialization of probabilities for each state for each observation
        initialize probability matrix
        first column of probability matrix:
            for each state: 
                initial probability of state k * emission probability for observation 0 in state k

        RECURSION: 
            For each remaining observation: 
                for each possible state:
                    for each possible previous state:
                        previous probability for preceding state  * transition probability * emission probability
                        add value to cumulative score counter (np.logaddexp)
                    store cumulative score for that cell in matrix, reset cumulative score
            return matrix

        TERMINATION:
            after all observations are iterated through
            sum values from last column --> overall probability

        return matrix and overall probability

    FUNCTION backward(self, observations):

        reverse order of observations
        call forward algorithm with reversed observations
        return matrix and overall probability

    FUNCTION forward_backward(self, observations):

        call forward algorithm
        call backward algorithm
        calculate average overall probability from forward and backward
        calculate marginal probability matrix for each state (fk + bk - overall_probability) in vector form
        return marginal probabilities
```

# Successes
There were a few key successes for this weeks project. We felt confident in our use of object oriented programming and we feel aptly set up to continue to build upon the objects we have created. Another highlight is the group work itself, we were able meet and discuss the algorithm conceptually, plan, and implement in an efficient manner. We ended up with a implementation we feel confident in moving forward with. 

For the second part, we came in with a good understanding of the forward, backward, and forward-backward algorithms. Our planning stage went well, as we were able to write detailed psuedocode for each function. This made the implementation very straightforward. We were also able to integrate our functions from this week into our HMM class quite easily and pull from the work we had done for viterbi. 

# Struggles
We initially struggled with the structure, or lack thereof, for this project. Setting up a notebook, and tackling this algorithm starting with nothing but our conceptual understanding required some teamwork and open discussion. We also had some trouble with the `_traceback()` function, keeping track of the correct column indices and what states to add to the best path. In paper, looking at the diagram in class, it seemed pretty straight-forward, but when we started coding it was a little more challenging. However, after thinking about it and osme trial-and-error, we were able to get it to work. 

For the second part of the project, using `np.logaddexp()` was a litte confusing at first. We weren't sure if using a cumulative score counter and repeatedly adding using the function would cause issues, and were struggled identifying when we were in log-space and when we weren't. However, after consulting the official documentation and some trial and error, we got it to work as expected. 

# Personal Reflections
## Group Leader
Spencer: I enjoyed this project, having no set structure and guidance as I think it emphasized the need for teamwork and made us leverage the skills we already have. I am looking forward to the following weeks working with Marcos and Sneha!

I found project 9 to be a bit easier to implement, but an apt follow up to viterbi. I appreciate our object oriented approach for project 8 as I think it set us up for success in this project. We ended up using initialization methods used in viterbi for our forward and backward algorithms, which made the impact of this approach clear!

## Other member
Sneha: Spencer and Marcos are both great partners and it was easy to meet and plan our implementation. We spent a good amount of time on the planning stage to make sure we understood each step of the algorithm and the scope of each of our functions. I think the most difficult concept for me was the traceback, specifically keeping track of the indices and mentally tracking that i and the actual matrix column are always one apart. It was also a bit intimidating at first not having a structured notebook and making sure that our code was implemented in such a way that we could add on to it in the next few weeks, but I think our group handled that well in planning. 

After understanding and implementing viterbi, I thought the second week of this project seemed pretty straightforward. One thing I struggled with conceptually at first was where the initial probabilities get added for the backwards algorithm, however after going through psuedocode, this became much more clear. I am proud of how we set up our project from the first week because we were able to intergrate our code for this week and use functions from last week pretty seamlessly. 

Marcos: I think our team worked very well given the lack of structure for this project. Just having a general idea of what we had to implement, with a few given data structures was a little daunting at first. However, during our planning meetings we were able to narrow down what the problem we had to implement was, and then how to actually do it. Having such a clear plan made the coding much easier. I really like how we decided to go with the class object approach, which allows us to have different methods and can continue to expand it with the coming weeks' projects. I think the hardest part was the `_traceback()` function, as I had a pretty solid idea of how it worked based on the table we saw in class, but it was challenging to translate it into actual code, keeping track of the correct indices. Despite that, we were able to overcome this and implemented the Viterbi algorithm correctly.

The second phase of the project was, in my opinion, slightly easier. Conceptually, the forward algorithm was intuitive and easy to implement. The only thing we really struggled with was understanding how to use `np.logaddexp()` correctly, as we got slightly confused about when we were in log-space and when we weren't. After consulting the official documentation and some trial and error, we got it to work as expected. 

# Generative AI Appendix
As per the syllabus
