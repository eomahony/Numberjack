#! /usr/bin/env python

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

    def same(category1, value1, category2, value2):
        return [(category1[i] == value1) == (category2[i] == value2) for i in range(Houses)]

    def next_to(category1, value1, category2, value2):
        constraints = [(category1[i] == value1) <=
                       ((category2[i-1] == value2) | (category2[i+1] == value2)) for i in range(1, Houses-1)]
        constraints.append((category1[0] == value1) <= (category2[1] == value2))
        constraints.append((category1[4] == value1) <= (category2[3] == value2))
        return constraints

    model += AllDiff(Nationalities)
    model += AllDiff(Colours)
    model += AllDiff(Animals)
    model += AllDiff(Drinks)
    model += AllDiff(Cigarettes)

    model += same(Nationalities, 'English', Colours, 'Red')
    model += same(Nationalities, 'Spanish', Animals, 'Dog')
    model += same(Drinks, 'Coffee', Colours, 'Green')
    model += same(Nationalities, 'Ukrainian', Drinks, 'Tea')
    model += same(Cigarettes, 'OldGold', Animals, 'Snail')
    model += same(Cigarettes, 'Kools', Colours, 'Yellow')
    model += same(Cigarettes, 'LuckyStrike', Drinks, 'OrangeJuice')
    model += same(Nationalities, 'Japanese', Cigarettes, 'Parliaments')

    model += Drinks[2] == 'Milk'
    model += Nationalities[0] == 'Norwegian'

    model += Colours[0] != 'Green'
    model += [(Colours[i] == 'Green') == (Colours[i - 1] == 'Ivory') for i in range(1, Houses)]

    model += next_to(Cigarettes, 'Chesterfields', Animals, 'Fox')
    model += next_to(Cigarettes, 'Kools', Animals, 'Horse')
    model += next_to(Nationalities, 'Norwegian', Colours, 'Blue')

    return Nationalities, Colours, Animals, Drinks, Cigarettes, Houses, model


def solve(param):
    Nationalities, Colours, Animals, Drinks, Cigarettes, Houses, model = get_model()

    solver = model.load(param['solver'])
    solver.solve()

    for i in range(Houses):
        print('The tenant of house', str(i+1), 'is', Nationalities[i], 'likes', Colours[i], 'has a', Animals[i], 'drinks', Drinks[i], 'and smokes', Cigarettes[i])


default = {'solver': 'Mistral'}


if __name__ == '__main__':
    param = input(default)
    solve(param)
