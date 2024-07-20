# Q1
class SuffixTrieNode:
    def __init__(self):
        """
            Function Description:
            Creates a new node to be used in SuffixTrie.

            Time and Space Complexity: O(1)
        """
        self.children = [None] * 4  # List to store references to child nodes. Size 4 to represent 'A,B,C,D'.
        self.indexes = [] # List to store indexes of substrings that this node makes up.

class SuffixTrie:
    def __init__(self, genome):
        """
            Function Description:
            Creates a suffix trie of all suffixes of 'genome', which stores the indexes of the starting position
            of these suffixes.

            Input: 
            genome, string variable.

            Time Complexity:
            O(N^2), where N is the number of characters in genome.
            Analysis:
            The outer loop 'i' iterates over the length of the genome, O(N).
            The inner loop 'j' iterates for (N+1)/2 times. 
            This is calculated from: i=0, 'j' iterates N times, i=1, 'j' iterates N-1 times, ... i=N-1, 'j' iterates 1 time.
            As a result, O(N * (N+1)/2) = O(N^2) when considering dominant elements.

            Space Complexity:
            Input and Auxiliary: O(N^2), N is the number of characters in genome.
            A new node is created for each unique substring. This is created in the inner loop, which runs for N*(N+1)/2 times
            (calculated previously in the time complexity). As a result, there are a maximum possible N*(N+1)/2 number
            of substrings, and thus number of nodes.
            The total number of substrings corresponds to the total number of indexes required to be stored. 
            This can be calculated from: i=0, N possible substrings starting at index 0, i=1, N-1 substrings starting at index 1
            ... i=N, 1 substring starting at index N.
            As a result, N*(N+1)/2 represents the space required for storing indexes (SuffixTrieNode.indexes) across all nodes.
            Thus, calculating complexity: O(2 * N*(N+1)/2) = O(2* N^2) = O(N^2) space required.
        """
        self.root = SuffixTrieNode()
        n = len(genome)
        
        for i in range(n): # Outer loop iterates over all the starting positions of suffixes. Runs for N times. 
            current_node = self.root
            # Inner loop iterates over each character in each suffix (this generates the unique substrings). 
            for j in range(i, n): # Runs for the length of the substring, N-i times.
                char_index = self.get_child_index(genome[j]) # Returns the numeric index to represent the current character.
                if not current_node.children[char_index]: # Creates a new node if child node does not exist.
                    current_node.children[char_index] = SuffixTrieNode()
                current_node = current_node.children[char_index]
                current_node.indexes.append(i) # Appends the starting index of the substring represented by 'i'.

    def search(self, string):
        """
            Function Description:
            Returns a list containing all the starting indexes of 'string' in self.genome.

            Input:
            string, string variable

            Output:
            list, containing all the start indexes of 'string'.

            Time Complexity:
            O(m), where m is the length of the string.
            Analysis:
            The worst case scenario, the loop will have to iterate through all characters in 'string'.
            All other operations are constant time, resulting in a time complexity of O(m), m is the length of 'string'.

            Space Complexity:
            Auxiliary: O(1) Input: O(I), I is the size of indexes.
            Analysis:
            No additional space is required by the function. 
            Input size of O(I) represents the size of current_node.indexes, which was already created in __init__.
        """
        current_node = self.root
        for char in string:
            char_index = self.get_child_index(char)
            if not current_node.children[char_index]: # The substring does not exist
                return []
            current_node = current_node.children[char_index]
        return current_node.indexes

    def get_child_index(self, char):
        """
            Function Description:
            Convert 'char' into a numeric index. A=0, B=1, C=2, D=3

            Input:
            char, string variable

            Output:
            integer, representing the index of the char.
        """
        return ord(char) - ord('A')


class OrfFinder:
    def __init__(self, genome):
        """
            Function Description:
            Creates a suffix trie containing all suffixes from 'genome'.

            Approach:
            This function creates a suffix trie to store all suffixes of 'genome'. This suffix trie contains nodes which
            each represent a character in a path of the trie. The path from the root to a node represents a substring of genome.
            Example suffix trie for genome 'ABC'.
          (root)
            ├── A (indexes: [0])          (substring 'A' begins at index 0)
            │   ├── B (indexes: [0])      (substring 'AB' begins at index 0)
            │   │   ├── C (indexes: [0])  (substring 'ABC' begins at index 0)
            ├── B (indexes: [1])          (substring 'B' begins at index 1)
            │   ├── C (indexes: [1])      (substring 'BC' begins at index 1)
            ├── C (indexes: [2])          (substring 'C' begins at index 2)
            Overall, this suffix trie allows for fast retrieval when searching for the starting indexes of a substring.

            Time Complexity:
            O(N^2), N is the length of the genome. 
            Analysis: See __init__ in class SuffixTrie.

            Space Complexity:
            Input and Auxiliary: O(N^2), N is the length of the genome. 
            Analysis: See __init__ in class SuffixTrie.
        """
        self.genome = genome
        self.suffix_trie = SuffixTrie(genome)

    def find(self, start, end):
        """
            Function Description:
            Returns a list containing all possible substrings of genome which begin with 'start' as a prefix and end with 'end' as a suffix.

            Approach Description:
            Searches the suffix tree initialised in __init__ to find the indexes of where 'start' and 'end' occur. 
            This returns an ordered list of the indexes of 'start' and 'end', from low to high. (Represented as start_indexes and end_indexes)
            Handles some edge cases, where 'end' is not found in the genome or the first occurrence of 'start' is after the last occurrence
            of 'end'.
            If the edge cases are not reached, then the function will iterate through the indexes of 'start' and 'end' which will produce
            a valid string and performs a list slice before appending it into outputs. 
            Invalid strings where the start index is at a position after the end index is avoided through checks across each iteration 
            for start_index and end_index where the next index for start and end is checked to see if a valid substring will be obtained.
            The start indexes iterate from min to max, whilst end indexes iterate from max to min to allow the remaining elements to be
            skipped entirely, avoiding redundant checks.

            Checking the next end index-
            If the next end index is found to not be at a position that is after the final character of the current start index,
            then the next substring will not be a valid as the 'end' prefix must be at a position after the 'start' prefix with no overlapping. 
            As a result, subsequent end indexes in the end_index list will also not produce valid substrings since we are iterating from 
            maximum to minimum, as they too will also be at an index less than the index of the current start index. 
            Thus, we break the loop for end indexes to avoid checking the redundant end indexes which will not produce a valid substring.

            Checking the next start index-
            If the next start index is at a position following the final (max) end index, then no substring can exist with the 'start'
            prefix of this index given that all positions of the 'end' suffix are found before it.
            As a result, subsequent start indexes in the start_index list will not produce valid substrings, since we are iterating 
            from minimum to maximum, as they too will also be at a index greater than the final index for 'end'.
            Thus, we break the loop for start indexes to avoid checking the redundant start indexes.

            These checks for next end index and next start index prevents the function from iterating over all the possible indexes in start_index
            and end_index, which would have resulted in a complexity of O(len(start_index) * len(end_index)).
            Instead, this approach ensures that only the indexes which will produce a valid string are iterated through, meeting the 
            complexity requirements.

            Input:
            start, string variable
            end, string variable

            Output:
            list, list of strings containing all the substrings in 'genome' that start with 'start' and end with 'end'

            Time Complexity:
            O(T + U + V)
            T representing the length of 'start', U representing the length of 'end', V representing the number of characters in the outputs list.
            Analysis:
            Searching the suffix tree for the string 'start' runs in O(T) time, searching for 'end' runs in O(U) time.
            In all cases, the nested loop will run for O(X) time, where X represents the number of substrings in the output list,
            as the loop will break if the next start index or end index will produce an invalid substring.
            In each iteration of the nested loop, an O(m) list slice is performed to create the valid substring which is appended to output O(1), 
            where m represents the length of that substring. 
            Since all substrings for where the list slice is performed will make up the final output list, the sum of all O(m) list slices
            across all iterations will add up to O(V), representing the number of characters in the final output list.
            As a result, summing these complexities will yield O(T + U + X + V). However, as X will always be less than V (the number of
            substrings in output will always be less than the number of characters in output), X will always be dominated by V.
            Thus, when considering the dominant elements, O(T + U + V) time complexity is reached.

            Space Complexity: 
            Input space complexity: O(N + V), where N is the length of genome. 
            Auxiliary space complexity: O(V), where V is the number of characters in the output list.
            Analysis:
            In the worst case scenario, the sum of the lists start_indexes and end_indexes would store N elements, where N is the total
            number of characters in genome. This would occur when the genome is only made up of 'start' or 'end', such as when 
            genome = 'AAABBB', and find('A', 'B') is called. In some cases, N > V such as when there are multiple indexes where 'start' 
            and 'end' are found at, but only a single valid substring in the output V. But in other cases, V > N, where there may be only
            1 index of 'start' and 'end' found, but the substring consists of multiple characters. As such, input space complexity can be
            represented as O(N + V) given that it is dependant on the inputs as to which is dominant. 
            The size of the lists start_indexes and end_indexes does not affect auxiliary space, as these do not require extra space, but
            are simply references to the already existing lists which were created within __init__. As a result, the auxiliary space for
            start_indexes and end_indexes is O(1). Outputs represents the extra auxiliary space required by the function, which is always
            O(V), where V is the number of characters found in the output list, resulting in O(V) auxiliary space.
        """
        start_indexes = self.suffix_trie.search(start) # Ordered list of indexes where 'start' occurs in the genome from min to max. O(T) time
        end_indexes = self.suffix_trie.search(end) # Ordered list of indexes where 'end' occurs in the genome from min to max. O(U) time

        outputs = []
        
        # Case when end is not found in the genome.
        # Returns outputs (empty list) to avoid iterating over the possible start and end indexes as no valid substring can exist.
        if len(end_indexes) == 0: # O(1) time
            return outputs

        # Case when the first occurrence of 'start' is after the last occurrence of 'end'.
        # This means that no substring exists that begins with 'start' and ends with 'end'.
        # Returns outputs (empty list) to avoid iterating over the possible start and end indexes as no valid substring can exist.
        if end_indexes[-1] <= start_indexes[0] + (len(start)-1): # -1 represents the final index in end_indexes, which is the greatest index. O(1) time.
            return outputs
        
        # Iterating through start_indexes from beginning to end (min to max), with the nested loop 
        # iterating through end_indexes in reverse (max to min) to ensure only valid indexes are checked.
        # In total, x number of iterations will occur, where x is the number of valid substrings in the output 
        # as it will break the loop if the next start or end index will result in an invalid substring.
        for i in range(len(start_indexes)):
            for j in range(1, len(end_indexes) + 1):
                # List slices the valid substring beginning with 'start' and ending with the final character of 'end'.
                outputs.append(self.genome[start_indexes[i]:end_indexes[-j] + len(end)]) # O(m), m is the number of characters in the list slice.
                if j < len(end_indexes): # If not the last element
                    # If the next end index is before or equal to the end of the start index, 
                    # break the loop as the remaining end indexes won't be valid. 
                    if end_indexes[-j-1] <= start_indexes[i] + (len(start) - 1): # O(1)
                        break
            if i < len(start_indexes) - 1: # If not the last element
                # If the maximum end index is before or equal to the end of the next start index, 
                # break the loop to avoid iterating remaining start indexes as they won't be valid.
                if end_indexes[-1] <= start_indexes[i + 1] + (len(start) - 1): # O(1) 
                    break

        return outputs