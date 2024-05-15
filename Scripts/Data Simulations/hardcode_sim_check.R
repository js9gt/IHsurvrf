
### simulation check


#### 1: check rates: failure, treatment, censoring

# failure coefficient
a1 = -6
b1 = -0.3
c1 = -0.5
z1 = -0.025
p1 = 0.02
g1 = -0.02
h1 = 0.008
r1 = 0.02

state1 = 7.770589e-01
state2 = 5.010484e-01
# ratio of input for max stages
nstages = 0.3
cumulative.length = 0.07
action = 1
prior.visit.length = 50

rate_failure = exp(
  a1 + b1 * state1 + c1*state2 + z1 * nstages + p1 * cumulative.length + g1 * action + h1*prior.visit.length +
    r1 * action * state1*state2 * nstages * cumulative.length*prior.visit.length
)

set.seed(127)
rexp(n = 1, rate = rate_failure)


## time to next visit coefficient


a2 = -3
b2 = -0.3
z2 = -0.015
p2 = 0.025
g2 = -0.2
h2 = 0.008
r2 = 0.01

## censoring coefficients
a3 = -8
b3 = -0.3
z3 = -0.04
p3 = 0.02
g3 = -0.2
h3 = 0.008
r3 = 0.01

#### 2: check new state

## this value should be between 0,1 since we are drawing from uniform distribution
next_state = 0.075440487
prev_state =1.35378074
rho = 0.75
# action-specific effects 
D0 = 0
D1 = 1
g = 0.25
action = 1

denom = ((1-(rho)^2)^0.5)*0.5*g


  
(( (next_state - rho*prev_state - ifelse(prev_state > 0.5, 1, 0)*(D0 - action*D1) - 
  ifelse(prev_state < 0.5, 1, 0)*(D0^2 + action*(D1)^2)))/ denom ) 




