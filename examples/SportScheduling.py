from Numberjack import *

def model_round_robin_tournament(nbTeams):
    print nbTeams
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

def solve_round_robin_tournament(param):
    (team, Weeks, Periods, model) = model_round_robin_tournament(param['teams'])
    solver = model.load(param['solver'])
    if solver.solve():
        print_schedule(team, Weeks, Periods)
    else:
        print 'no such tournament'
    print 'Nodes:', solver.getNodes(), ' Time:', solver.getTime()

def print_schedule(team, Weeks, Periods):
    print '           ', ''.join([('week '+str(w+1)).ljust(9) for w in Weeks])
    for p in Periods:
        print ('period '+str(p)+':').ljust(10), ''.join([(str(team[w][p]['home']).rjust(2)+' vs '+str(team[w][p]['away']).ljust(2)).ljust(9) for w in Weeks]) 

solve_round_robin_tournament(input({'solver':'Mistral', 'teams':8}))




