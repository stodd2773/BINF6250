import unittest
from project01 import read_file, parse_line

class TestProject01(unittest.TestCase): 

    def test_DiseaseRare(self):
        '''
        Tests if disease is counted in dictionary when disease is rare.
        Expecting dictionary with only rare diseases.
        '''

        expectation = {'rare_disease': 1}
        result = read_file('raredisease.vcf')
        self.assertEqual(expectation, result)
    
    def test_DiseaseNotRare(self):
        '''
        Testing against not rare disease types 
        '''
        expectation = {}
        result = read_file('notraredisease.vcf')
        self.assertEqual(expectation, result)

    def test_MultDiseaseRare(self):

        '''
        Testing if dictionary includes multiple diseases associated with single rare allele 
        '''

        expectation = {'rare_disease_1': 1, 'rare_disease_2': 1}
        result = read_file('multraredisease.vcf')
        self.assertEqual(expectation, result)

    def test_NotProvided_MultDiseaseRare(self):

        '''
        Testing if dictionary includes second disease associated with single rare allele, if the first disease is 'not provided' 
        '''

        expectation = {'rare_disease_2': 1}
        result = read_file('notprovided_multraredisease.vcf')
        self.assertEqual(expectation, result)
        
        
unittest.main()        

