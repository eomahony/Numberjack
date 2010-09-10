from Numberjack import *

def get_model(nbTeams):
    Teams   = range(nbTeams)
    Weeks   = range(nbTeams-1)
    EWeeks  = range(nbTeams)
    Periods = range(nbTeams/2)
    Games   = range(nbTeams*nbTeams)
    Slots   = ('home', 'away')

    occur = dict([(i,(2,2)) for i in range(nbTeams)])

    team = [[dict([(s,Variable(Teams)) for s in Slots]) for p in Periods] for w in EWeeks]
    game = [[Variable(Games) for p in Periods] for w in Weeks]

    model = Model(
        AllDiff( [game[w][p] for p in Periods for w in Weeks] ), # each game is played once
        [AllDiff( [team[w][p][s] for s in Slots for p in Periods] ) for w in EWeeks], # each team play at most once a week
        [Gcc( [team[w][p][s] for s in Slots for w in EWeeks], occur ) for p in Periods ], # each team play twice per period

        [team[w][p]['home']*nbTeams + team[w][p]['away'] == game[w][p] for w in Weeks for p in Periods] #channelling teams & games
        )

    #print model

    return (team, Weeks, Periods, model)

def solve(param):
    (team, Weeks, Periods, model) = get_model(param['teams'])
    solver = model.load(param['solver'])
    solver.setVerbosity(param['verbose'])
    solver.setTimeLimit(param['tcutoff'])
    solver.solve()

    out = ''
    if solver.is_sat():
        out = str(print_schedule(team, Weeks, Periods))
    out += ('\nNodes: ' + str(solver.getNodes()))
    return out


def print_schedule(team, Weeks, Periods):
    out = '            ' + ''.join([('week '+str(w+1)).ljust(9) for w in Weeks]) + '\n'
    for p in Periods:
        out += ('period '+str(p)+':').ljust(10) + ' ' + ''.join([(str(team[w][p]['home']).rjust(2)+' vs '+str(team[w][p]['away']).ljust(2)).ljust(9) for w in Weeks]) + '\n'
    return out


solvers = ['Mistral', 'MiniSat', 'SCIP', 'Walksat']
default = {'solver':'Mistral', 'teams':8, 'verbose':1, 'tcutoff':3}

if __name__ == '__main__':
    param = input(default) 
    print solve(param)



