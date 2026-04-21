import numpy as np            # Workhorse
import numpy.random as nr     # Setting up random distributions
import pandas as pd           # Simple convergence checks
from itertools import product # For robust iteration
from copy import deepcopy     # For convergence checking
from array import array       # Faster base Python iteration

try:
    # Special Json library that isn't only more efficient, 
    # but also allows for setting floating point precision on print
    import ujson as json

    def to_json(data, precision = 2):
        """Wrapper function to output dict to json with specific floating-point
        precision.

        Args:
            data (dict): some HMM probability matrix
            precision (int): Setting for floating-point precision (default: 2)
        
        Returns:
            (json): JSON representation of the dictionary
        """
        return json.dumps(data, indent = 4, double_precision = precision)

except ImportError:
    import json

    def to_json(data, precision = 2):
        """Wrapper function to output dict to json with specific floating-point
        precision.

        Note: Since the base `json.dumps` doesn't provide double_precision,
        this goes through an extra encoder step to handle the floats

        Args:
            data (dict): some HMM probability matrix
            precision (int): Setting for floating-point precision (default: 2)

        Returns:
            (json): JSON representation of the dictionary
        """
        out_json = json.dumps(data)
        loaded_json = json.loads(out_json, parse_float=lambda o: f'{float(o):.2g}')
        return json.dumps(loaded_json, sort_keys = True, indent = 4)

# Gil-ify
np.set_printoptions(precision = 2)
pd.set_option('display.precision', 2)

class BaseHMM(object):
    """Base class for HMM objects

    Class for holding HMM parameters and to allow for implementation of
    functions associated with HMMs

    Attributes:
        alphabet (str): The emissions used in the HMM (default: 'ACGT')
        hidden_states (str): The hidden states within the HMM (default: None)
        init_probs (dict of floats): β probabilities for initial steps (default: None)
        trans_probs (dict of dict of floats): Transition probabilities from one state to another given a state (default: None)
        emit_probs (dict of dict of floats): Emission probabilities of a letter given a state (default: None)
    """

    __all__ = ['alphabet', 'hidden_states', 'emit_probs', 'trans_probs', 'init_probs']

    def __init__(self, alphabet = "ACGT", hidden_states = None, init_probs = None, 
                 trans_probs = None, emit_probs = None, seed = None, precision = 2, tolerance = 1e-10):
        """Instantiate the HMM

        Args:
            alphabet (str): The emissions used in the HMM (default: 'ACGT')
            hidden_states (str): The hidden states within the HMM (default: None)
            init_probs (dict of floats): β probabilities for initial steps (default: None)
            trans_probs (dict of dict of floats): Transition probabilities from one state to another given a state (default: None)
            emit_probs (dict of dict of floats): Emission probabilities of a letter given a state (default: None)
            seed (int): To set the random seed for numpy (default: None)
            precision (int): Floating-point precision (default: 2)
            tolerance (number): The acceptable tolerance between baum-welch iteration to determine if convergence has occured (default: 1e-10)
        """

        # Floating point precision
        self._precision = precision

        # For convergence checking
        self._tolerance = tolerance

        # Need at least the hidden states to initialize
        if hidden_states is None:
            raise ValueError('Hidden states must be provided')
        self.hidden_states = hidden_states
        self.alphabet = alphabet

        # Cursory initialization that allows users to set some prior
        self.init_probs = init_probs
        self.trans_probs = trans_probs
        self.emit_probs = emit_probs
        
        # If the user doesn't provide any probabilities, randomize some stuff
        self._initialize_random(seed)
    
    def _initialize_random(self, seed):
        """[PRIVATE] Used to completely initialize the HMM if any probability matrices
        provided. Will not overwrite given probabilities.

        Args:
            seed (int): To set the random seed for numpy
        """
        nr.seed(seed)
        if self.init_probs is None:

            # These aliases help keep line width down during the iteration and allows 
            # for simple dictionary comprehension
            states = self.hidden_states
            i_probs = nr.dirichlet(np.ones(len(self.hidden_states)))

            # Make a dictionary of floats (that sum up to 1) for all states in the HMM
            init = {state: np.float64(i_prob) for state, i_prob in zip(states, i_probs)}
            self.init_probs = init

        if self.trans_probs is None:

            # These aliases help keep line width down during the iteration and allows 
            # for simple dictionary comprehension
            states = product(self.hidden_states, self.hidden_states)
            t_probs = np.nditer(nr.dirichlet(np.ones(len(self.hidden_states)), size=len(self.hidden_states)))

            # Make a dictionary for each state that contains a dictionary of floats
            # that sum to 1 and represent each of the states used in the HMM
            trans = {}
            for (state, next_state), t_prob in zip(states, t_probs):
                trans.setdefault(state, {}).update({next_state: np.float64(t_prob)})
            self.trans_probs = trans

        if self.emit_probs is None:

            # These aliases help keep line width down during the iteration and allows 
            # for simple dictionary comprehension
            states_letters = product(self.hidden_states, self.alphabet)
            e_probs = np.nditer(nr.dirichlet(np.ones(len(self.alphabet)), size=len(self.hidden_states)))

            # Make a dictionary for each state that contains a dictionary of floats
            # that sum to 1 and represent each of the emissions used in the HMM
            emit = {}
            for (state, letter), e_prob in zip(states_letters, e_probs):
                emit.setdefault(state, {}).update({letter: np.float64(e_prob)})
            self.emit_probs = emit

    # All the `@property` and `@''.setter` methods are doing are allowing
    # for isolation and encapsulation of the internal datasets

    @property
    def init_probs(self):
        return self._initial

    @init_probs.setter
    def init_probs(self, init_probs):
        if init_probs is None:
            self._initial = None
        elif isinstance(init_probs, dict):
            if len(init_probs) != len(self.hidden_states):
                raise ValueError('Initial probabilites must be the length of the number of hidden states')
            if not np.isclose(sum(init_probs.values()), 1):
                raise ValueError('Initial probabilites must sum to 1')
            self._initial = init_probs
        else:
            raise SyntaxError('Initial probabilities must be None or a dictionary')

    @property
    def trans_probs(self):
        return self._transition

    @trans_probs.setter
    def trans_probs(self, trans_probs):
        if trans_probs is None:
            self._transition = None
        elif isinstance(trans_probs, dict):
            if len(trans_probs) != len(self.hidden_states):
                raise ValueError('Transition probabilites must be a square matrix')
            if len(trans_probs[self.hidden_states[0]]) != len(self.hidden_states):
                raise ValueError('Transition probabilites must be a square matrix')
            if not np.allclose([sum(trans_probs[state].values()) for state in self.hidden_states], 1):
                raise ValueError('Transition probabilites must sum to 1 along a given axis')
            self._transition = trans_probs
        else:
            raise SyntaxError('Transition probabilities must be None or a dictionary')

    @property
    def emit_probs(self):
        return self._emission

    @emit_probs.setter
    def emit_probs(self, emit_probs):
        if emit_probs is None:
            self._emission = None
        elif isinstance(emit_probs, dict):
            if len(emit_probs) != len(self.hidden_states):
                raise ValueError('Emission probabilites must be length of hidden states by length of alphabet')
            if len(emit_probs[self.hidden_states[0]]) != len(self.alphabet):
                raise ValueError('Emission probabilites must be length of hidden states by length of alphabet')
            emit = pd.DataFrame.from_dict(emit_probs).T
            emit.columns = list(self.alphabet)
            if not np.allclose([sum(emit_probs[state].values()) for state in self.hidden_states], 1):
                raise ValueError('Emission probabilites must sum to 1 along a given axis')
            self._emission = emit_probs
        else:
            raise SyntaxError('Emission probabilities must be None or a dictionary')

    @property
    def hidden_states(self):
        return self._hidden_states

    @hidden_states.setter
    def hidden_states(self, hidden_states):
        if isinstance(hidden_states, str):
            self._hidden_states = hidden_states
        elif isinstance(hidden_states, (tuple, list)):
            self._hidden_states = ''.join(hidden_states)

    @property
    def alphabet(self):
        return self._alph

    @alphabet.setter
    def alphabet(self, alphabet):
        if isinstance(alphabet, str):
            self._alph = alphabet
        elif isinstance(alphabet, (tuple, list)):
            self._alph = ''.join(alphabet)

    def __str__(self):
        out_text = [f'Alphabet: {self.alphabet}',
                    f'Hidden States: {self.hidden_states}',
                    f'Initial Probabilities: {to_json(self.init_probs, self._precision)}',
                    f'Transition Probabilities: {to_json(self.trans_probs, self._precision)}',
                    f'Emission Probabilities: {to_json(self.emit_probs, self._precision)}']
        return '\n'.join(out_text)

    @classmethod
    def __dir__(cls):
        return cls.__all__

    def __eq__(self, other):
        """Solely used by the `baum_welch` function to check for convergence"""
        if np.allclose(pd.Series(self.init_probs), pd.Series(other.init_probs), atol=self._tolerance):
            if np.allclose(pd.DataFrame(self.trans_probs), pd.DataFrame(other.trans_probs), atol=self._tolerance):
                if np.allclose(pd.DataFrame(self.emit_probs), pd.DataFrame(other.emit_probs), atol=self._tolerance):
                    return True
        return False


class HMM(BaseHMM):
    """Main class for HMM objects

    Class for holding HMM parameters and to allow for implementation of
    functions associated with HMMs

    Attributes:
        alphabet (str): The emissions used in the HMM (default: 'ACGT')
        hidden_states (str): The hidden states within the HMM (default: None)
        init_probs (dict of floats): β probabilities for initial steps (default: None)
        trans_probs (dict of dict of floats): Transition probabilities from one state to another given a state (default: None)
        emit_probs (dict of dict of floats): Emission probabilities of a letter given a state (default: None)
    """

    __all__ = ['alphabet', 'hidden_states', 'emit_probs', 'trans_probs', 'init_probs',
               'viterbi','forward', 'backward', 'forward_backward', 'baum_welch']

    @classmethod
    def __dir__(cls):
        return cls.__all__

    def forward(self, sequence):
        """ The forward algorithm for calculating probability of sequence given HMM

        Args:
            sequence (str): a valid set of emissions from the HMM

        Returns:
            forward_prob (float): probability of the the sequence using the forward algorithm
            forward (dict of array of float): the forward matrix of the probabilities of given bases at given positions and given states
        """
        # Initialize the lattice. I am using an array here because it is more
        # performant in iteration and has a smaller footprint than mutliple dictionaries
        forward = {state: array('d', [0] * len(sequence)) for state in self.hidden_states}

        # Since the forward algorithm starts at the beginning, use the first emission
        # and initial state probabilities to fill in the first column
        for state in self.hidden_states:
            forward[state][0] = self.init_probs[state] * self.emit_probs[state][sequence[0]]

        # This is where things get interesting: by taking the cartesian product of the
        # index position (along the sequence, starting from 1 since 0 is already filled in) and
        # The hidden states, we can condense a nested for-loop into a single line.
        for seq_idx, next_state in product(range(1, len(sequence)), self.hidden_states):

            # Since the forward algorithm just takes the sum across states from a given state, we need 
            # to tease each of the states out.
            # Furthermore, this is always using the data that has already been calculated from the left.
            for curr_state in self.hidden_states:
                forward[next_state][seq_idx] += forward[curr_state][seq_idx - 1] * self.trans_probs[curr_state][next_state]

            # Now that we have our positional sum, we use the emission probability for that state to update
            forward[next_state][seq_idx] = (self.emit_probs[next_state][sequence[seq_idx]] 
                                          * forward[next_state][seq_idx])

        # When all is done, the final probability of the sequence (based on the forward algorithm),
        # is the sum across both states at the end
        forward_prob = sum(forward[state][-1] for state in self.hidden_states)
        return forward_prob, forward 

    def backward(self, sequence):
        """ The backward algorithm for calculating probability of sequence given HMM

        Args:
            sequence (str): a valid set of emissions from the HMM

        Returns:
            backward_prob (float): probability of the the sequence using the backward algorithm
            backward (dict of array of float): the backward matrix of the probabilities of given bases at given positions and given states
        """
        # Initialize the lattice. I am using an array here because it is more
        # performant in iteration and has a smaller footprint than mutliple dictionaries
        backward = {state: array('d', [0] * len(sequence)) for state in self.hidden_states}

        # The backward algorithm starts at the end. Make that 1 and work from there
        for state in self.hidden_states:
            backward[state][-1] = 1

        # Like the forward algorithm: by taking the cartesian product of the
        # index position (along the sequence starting from the end since) and
        # The hidden states, we can condense a nested for-loop into a single line.
        rev_seq = range(len(sequence) - 1, 0, -1)
        for seq_idx, last_state in product(rev_seq, self.hidden_states):

            # Since the forward algorithm just takes the sum across states from a given state, we need 
            # to tease each of the states out
            # Furthermore, this is always using the data that has already been calculated from the right.
            for curr_state in self.hidden_states:
                backward[last_state][seq_idx-1] += backward[curr_state][seq_idx] * self.trans_probs[last_state][curr_state] * self.emit_probs[curr_state][sequence[seq_idx]]

        # When all is done, the final probability of the sequence (based on the backward algorithm),
        # is the sum across both states at the start (relative to left-to-right)
        backward_prob = sum(backward[state][0]*self.init_probs[state] * self.emit_probs[state][sequence[0]] for state in self.hidden_states)
        return backward_prob, backward

    def forward_backward(self, sequence):
        """ The forward-backward algorithm for calculating marginal posteriors given HMM

        Args:
            sequence (list): a list of valid emissions from the HMM

        Returns:
            posterior (list of dicts): all posteriors as a list
        """
        #Calculate forward and backward matrices
        Pf, f_matrix = self.forward(sequence)
        Pb, b_matrix = self.backward(sequence)

        # Generally speaking, most will only use either the forward or the backward
        # probability of the sequence...not both. However, these two probabilities
        # should be nearly the identical and I like being conservative. So I use the
        # average and everybody is included.
        P = (Pf + Pb)/2

        # Initialize the lattice. I am using an array here because it is more
        # performant in iteration and has a smaller footprint than mutliple dictionaries
        posterior = {state: array('d', [0] * len(sequence)) for state in self.hidden_states}

        # By using the cartesian product of the range of the sequence length and
        # the hidden states, I can condense a nested for-loop to a singe line
        for i, state in product(range(len(sequence)), self.hidden_states):
            posterior[state][i] = f_matrix[state][i] * b_matrix[state][i] / P  
        return posterior

    def viterbi(self, sequence):
        """ The viterbi algorithm for decoding a string using a HMM

        Args:
            sequence (str): a list of valid emissions from the HMM

        Returns:
            result (str): optimal path through HMM given the model parameters
                           using the Viterbi algorithm
        """
        def update_probs(base, previous):
            """Nested function used to keep track of the current probabilities and update the next

            Args:
                base (str): the current emission
                previous (dict of float): previous position's probabilities

            Returns:
                next_prob (dict of float): Next position's probabilities
                tb (dict of str): The traceback from current to previous origin
            """
            curr_prob = {} # Will caclculate our current position's probabilities
            next_prob = {} # Will contain the new position's probabilities
            tb = {}        # Contains the computed traceback as {current_state: previous} entries

            for next_state in self.hidden_states:
                for curr_state in self.hidden_states:
                    curr_prob[curr_state] = previous[curr_state] + np.log10(self.trans_probs[curr_state][next_state])

                # This max function acts as a argmax. This is because the key parameter can take a function
                # to determine how max is computed. Here we are telling it to base it on the values of the keys
                # and not the keys themselves and then return the key matching the max value
                origin = max(curr_prob, key=curr_prob.get)

                # We use that origin the next states probability based on the current emission
                next_prob[next_state] = np.log10(self.emit_probs[next_state][base]) + curr_prob[origin]
                tb[next_state] = origin
            return next_prob, tb

        def get_traceback(traceback, last_origin):
            """Nested function that parses the traceback dict and constructs the most optimal 
            path of states given the sequence

            Args:
                traceback (dict of str): The traceback of all positions to their origins
                last_origin (str): the max state from last position of the probability matrix

            Returns:
                tb (str):  the rest of the path, starting from last_origin
            """
            tb = ''

            # Reverse the traceback so that we start at the end
            for pos in reversed(traceback):
                # We already determine last_origin based on the final outcome
                # of the probability matrix
                prev_origin = pos[last_origin]

                # Keep adding to our sequence of optimal origins
                tb += prev_origin

                # Update for the next iteration
                last_origin = prev_origin
            return tb

        traceback = []

        first_base = sequence[0]

        # Start off by using the initial conditions and the first emission of the sequence
        previous = {state: np.log10(self.init_probs[state]) + np.log10(self.emit_probs[state][first_base]) for state in self.hidden_states}

        # Go through all other positions and keep track of the running total of probabilities
        for base in sequence[1:]:
            update_previous, update_tb = update_probs(base, previous)
            previous = update_previous
            traceback.append(update_tb)

        # Find the max state at the final position
        result = max(previous, key=previous.get)

        result += get_traceback(traceback, result)        

        # Since Traceback starts at the end and works forward, 
        # we need to reverse the result
        return result[::-1]

    def baum_welch(self, sequences, pseudocount = 1e-100):
        """Baum-Welch is an EM-algorithm that finds the maximum likelihood estimate of 
        the parameters of a HMM given a set of observed emission sequences.

        Note: Used when the user doesn't know all/any of the HMM's probabilities.

        Args:
            sequences (list of str): all the sequences used for training the HMM
            pseudocount (number): some pseudocount to prevent ZeroDivisionError (default: 1e-100)
        """
        def init_bw(pseudocount):
            """Initializes the pseudocount-filled probability matrices for init, trans, and emit

            Args:
                pseudocount (number): some pseudocount to prevent ZeroDivisionError

            Returns:
                init (dict of floats): Pseudocount-filled matrix for initial steps
                trans (dict of dict of floats): Pseudocount-filled matrix for transition probabilities from one state to another given a state
                emit (dict of dict of floats): Pseudocount-filled matrix for emission probabilities of a letter given a state
            """
            init = {state: pseudocount for state in self.hidden_states}
            trans = {state: {next_state: pseudocount for next_state in self.hidden_states} for state in self.hidden_states}
            emit = {state: {letter: pseudocount for letter in self.alphabet} for state in self.hidden_states}
            return init, trans, emit
        
        def proc_seq(outer, seq, init, trans, emit):
            """Processes a given sequence such that init, trans, and emit 
            are updated as specific emission are computed.

            Args:
                outer (float): the sum of observed sequence probabilities
                init (dict of floats): β probabilities for initial steps
                trans (dict of dict of floats): Transition probabilities from one state to another given a state
                emit (dict of dict of floats): Emission probabilities of a letter given a state

            Returns:
                outer (float): incremented sum of observed sequence probabilities
                init (dict of floats): Scaled β probabilities for initial steps
                trans (dict of dict of floats): Scaled Transition probabilities from one state to another given a state
                emit (dict of dict of floats): Scaled Emission probabilities of a letter given a state
            """
            prob_forward, forward = self.forward(seq)
            prob_backward, backward = self.backward(seq)

            # Fun Durbin step because, ummm...stats?
            # Generally speaking, most will only use either the forward or the backward
            # probability of the sequence...not both. However, these two probabilities
            # should be nearly the identical and I like being conservative. So I use the
            # average and everybody is included.
            outer += (prob_forward + prob_backward)/2

            # As I go through the sequence one emission at a time...
            for i, emission in enumerate(seq):
                # and visit each possible state that emission can be from...
                for state in self.hidden_states:
                    # I need to take into account the very first observation coming from the initial step
                    if i == 0:
                        init[state] += forward[state][i] * backward[state][i]
                    # and update the emission matrix for every emission observed at the given step
                    emit[state][emission] += forward[state][i] * backward[state][i]

                    if i == len(seq) - 1:
                        break # I have reached the end of the sequence wrt to transitions, so I stop

                    # Transitions are fun because it is always based on where it can go from where it is.
                    # Therefore we are always looking ahead
                    trans_emission = seq[i+1]
                    for next_state in self.hidden_states:
                        # Based on the 'future' state given our 'current' state, we use the probability of
                        # the next emission at the next position at the next state to update our transition
                        # matrix
                        trans[state][next_state] += (
                            forward[state][i] 
                            * self.trans_probs[state][next_state]
                            * self.emit_probs[next_state][trans_emission] 
                            * backward[next_state][i+1]
                        )

            return outer, init, trans, emit

        def scale_step(outer, init, trans, emit):
            """Scales all the probability matrices based on the sum of
            observed sequence probabilities.

            Args:
                outer (float): the sum of observed sequence probabilities
                init (dict of floats): β probabilities for initial steps
                trans (dict of dict of floats): Transition probabilities from one state to another given a state
                emit (dict of dict of floats): Emission probabilities of a letter given a state

            Returns:
                init (dict of floats): Scaled β probabilities for initial steps
                trans (dict of dict of floats): Scaled Transition probabilities from one state to another given a state
                emit (dict of dict of floats): Scaled Emission probabilities of a letter given a state
            """
            # Use that fun outer denominator we have been keeping track of. However,
            # if only one sequence is provided, the user could just use the posterior
            # (or forward_backward algorithm) instead.
            # For each of the matrices, we are essentially dividing all of our observed
            # probabilities by the probability of the sequences. If we don't we will greatly
            # underestimate our probabilities and likely hit an underflow issue
            for state in self.hidden_states:
                init[state] /= outer
                for letter in self.alphabet:
                    emit[state][letter] /= outer
                for other_state in self.hidden_states:
                    trans[state][other_state] /= outer

            # I think this is the Maximization step. However, since all of our probabilities
            # should just sum to 1 for any given state, this shouldn't have a large impact
            init_sum = sum(init.values())
            emit_sum = {state: sum(emit[state].values()) for state in self.hidden_states}
            trans_sum = {state: sum(trans[state].values()) for state in self.hidden_states}            

            for state in self.hidden_states:
                init[state] /= init_sum
                for letter in self.alphabet:
                    emit[state][letter] /= emit_sum[state]
                for other_state in self.hidden_states:
                    trans[state][other_state] /= trans_sum[state]
            return init, trans, emit

        converged = False
        count = 0

        while not converged:
            init, trans, emit = init_bw(pseudocount)

            # This is only used if there are multiple sequences
            # otherwise, forward-backward would be okay
            outer = 0 

            for seq in sequences:

                # Send each sequence to processing
                outer, init, trans, emit = proc_seq(outer, seq, init, trans, emit)

            # Now scale the probability matrices based on sum of observed sequence probabilities
            init, trans, emit = scale_step(outer, init, trans, emit)

            # Used for convergence checking
            old = deepcopy(self)

            # Update the model
            self.init_probs = init
            self.emit_probs = emit
            self.trans_probs = trans

            count += 1
            # Utilize that special HMM.__eq__ method
            if self == old:
                converged = True
                print(f'Converged after {count} iterations')

class pHMM(HMM):

    def read_fasta(self, path):

        # path = path_to_fasta
        # Read fasta file by lines
        # Skip headers
        # Add aligned sequences to list/matrix
        # return and create self.MSA
        pass

    def assign_state(self):

        # for each column (position) in self.MSA
            # Assign state to column based on alignment
                # state = match if >50% aligned
                # state = insertion if <50% aligned

        # return assignment matrix/list
        pass

    def label_position(self):
        
        # for each sequence in MSA
            # for each position in sequence 
                # Label position with residue and state
                    # match or insertion if there is a residue present
                    # match if global assignment is match
                    # insertion if global assignment is insertion
                    # deletion if '-' or symbol denoting gap

        # return array/matrix (row = sequence, col = position) with labels
        pass

    def background_dist(self):

        # take all observed residues and calculate prob for each residue
            # make empty background_dist_dict
            # for each residue: 
                # background_dist = (count of residue + pseudocount) / (total num of residues observed + 20 * pseudocount)
                # bacground_dist_dict[residue] = background_dist
        # return background_dist_dict
        pass

    def compute_emissions(self):

        # Use matrix with labeled positions for each sequence (self.label_position())
        # for state_index, residue_index in enumerate(labeled_martrix)
            
        
            # if state assignment is a match, execute build_emissions_from_msa for that column using only aligned residues (skipping gaps)
            #     residues = [seq[residue_index] for seq in msa 
            #     if seq[residue_index] != '-']

            #     def build_emissions(self, residues, state_index, emissions):
                        
            #         # Count frequencies
            #         counts = {}
            #         for r in self.alphabet:
            #             counts[r] = residues.count(r)
            #         total = sum(counts.values())
                    
            #         # Compute emission probabilities for match state
            #         self.emissions[f'M{state_index}'] = {
            #             r: counts[r] + self.pseudocount / total + 20 * self.pseudocount 
            #             for r in self.alphabet
            #         }
                    
            #         # Insertion state uses background distribution
            #         emissions[f'I{state_index}'] = self.background_dist()
                
            #         return self.emissions

            #     if state assignment is not match, skip
            pass
    
    def calc_trans_prob(self):

        # initialize with transition_counts with valid transitions: 
            # M_i -> M_i+1
            # M_i -> I_i 
            # M_i -> D_i+1
            # ...

        # for each labeled state in labeled_state_matrix.shape()
            # curr_state <- current labeled state
            # next_state <- next state in labeled state matrix

            # if curr_state to next_state transition is in transition count dictionary: 
                # increment counter for that transition by 1

        
        # for each key in the dictionary (original labeled state): 
            # sum all counts to get total
            # for each key (next state): 
                # value (transition count) + psuedocount / total counts + 20 * psuedocount
        # return self.trans_probs
        pass


    def evaluate(self, sequences):

        # for each seq in sequences: 
            # prob = pHMM.forward() to get overall prob
            # path = pHMM.viterbi() to get optimal state path for M, D, I
        # print prob and path variables 

        pass  