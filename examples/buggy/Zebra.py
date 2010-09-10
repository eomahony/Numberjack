#! /usr/bin/env python

import sys

from Numberjack import *

#-----Zebra Model--------------------

def get_model():
    model = Model()

    Nationalities = VarArray(5, ['English', 'Spanish', 'Ukrainian', 'Norwegian', 'Japanese'])
    Colours = VarArray(5, ['Red', 'Green', 'Ivory', 'Yellow', 'Blue'])
    Animals = VarArray(5, ['Dog', 'Fox', 'Horse', 'Zebra', 'Snail'])
    Drinks = VarArray(5, ['Coffee', 'Tea', 'Milk', 'OrangeJuice', 'Water'])
    Cigarettes = VarArray(5, ['OldGold', 'Kools', 'Chesterfields', 'LuckyStrike', 'Parliaments'])
    Houses = 5

    model.add( AllDiff(Nationalities) )
    model.add( AllDiff(Colours) )
    model.add( AllDiff(Animals) )
    model.add( AllDiff(Drinks) )
    model.add( AllDiff(Cigarettes) )


    def same( category1, value1, category2, value2 ):
        return [(category1[i] == value1) == (category2[i] == value2) for i in range(Houses)]

    model.add( same(Nationalities, 'English', Colours, 'Red') )
    model.add( same(Nationalities, 'Spanish', Animals, 'Dog') )
    model.add( same(Drinks, 'Coffee', Colours, 'Green') )
    model.add( same(Nationalities, 'Ukrainian', Drinks, 'Tea') )
    model.add( same(Cigarettes, 'OldGold', Animals, 'Snail') )
    model.add( same(Cigarettes, 'Kools', Colours, 'Yellow') )
    model.add( same(Cigarettes, 'LuckyStrike', Drinks, 'OrangeJuice') )
    model.add( same(Nationalities, 'Japanese', Cigarettes, 'Parliaments') )


    model.add( (Drinks[2] == 'Milk') )
    model.add( (Nationalities[0] == 'Norwegian') )

    model.add( (Colours[0] != 'Green') )
    model.add( [(Colours[i] == 'Green') == (Colours[i-1] == 'Ivory') for i in range(1,Houses)] )


    def next_to( category1, value1, category2, value2 ):
        constraints = [(category1[i] == value1) <= 
                       ((category2[i-1] == value2) | (category2[i+1] == value2)) for i in range(1, Houses-1)]
        constraints.append( (category1[0] == value1) <= (category2[1] == value2) )
        constraints.append( (category1[4] == value1) <= (category2[3] == value2) )
        return constraints

    model.add( next_to(Cigarettes, 'Chesterfields', Animals, 'Fox') )
    model.add( next_to(Cigarettes, 'Kools', Animals, 'Horse') )
    model.add( next_to(Nationalities, 'Norwegian', Colours, 'Blue') )

    return Nationalities,Colours,Animals,Drinks,Cigarettes,model


def solve(param):

    Nationalities,Colours,Animals,Drinks,Cigarettes,model = get_model()

    print Nationalities
    print Colours

    solver = model.load( param['solver'] )

    solver.solve()

    print Nationalities
    print Colours

    for i in range(Houses):
        print 'The tenant of house', str(i+1), 'is', Nationalities[i], 'likes', Colours[i], 'has a', Animals[i], 'drinks', Drinks[i], 'and smokes', Cigarettes[i]


solvers = ['Mistral', 'MiniSat', 'SCIP', 'Walksat']
default = {'solver':'Mistral'}

if __name__ == '__main__':
    param = input(default) 
    solve(param)



