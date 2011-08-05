import unittest
import doctest
import cloudtracker

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(cloudtracker))
        return tests    
