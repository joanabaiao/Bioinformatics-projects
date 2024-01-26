# Viterbi Algorithm 

# five elements for HMM
states = ('Healthy', 'Fever')

observations = ('normal', 'cold', 'dizzy')

start_probability = {'Healthy': 0.6, 'Fever': 0.4}

transition_probability = {'Healthy': {'Healthy': 0.7, 'Fever': 0.3},
                          'Fever': {'Healthy': 0.4, 'Fever': 0.6}}

emission_probability = {'Healthy': {'normal': 0.5, 'cold': 0.4, 'dizzy': 0.1},
                        'Fever': {'normal': 0.1, 'cold': 0.3, 'dizzy': 0.6}}


def Viterbit(obs, states, s_pro, t_pro, e_pro):
    
    # init path: path[s] represents the path ends with s
    path = {s: [] for s in states}
    curr_pro = {}
    for s in states:
        curr_pro[s] = s_pro[s]*e_pro[s][obs[0]]

    for i in range(1, len(obs)):
        last_pro = curr_pro
        curr_pro = {}
        for curr_state in states:
            max_pro, last_sta = max(((last_pro[last_state]*t_pro[last_state][curr_state]
                                      * e_pro[curr_state][obs[i]], last_state) for last_state in states))
            curr_pro[curr_state] = max_pro
            path[curr_state].append(last_sta)

    # Final largest probability
    max_pro = -1
    max_path = None
    for s in states:
        path[s].append(s)
        if curr_pro[s] > max_pro:
            max_path = path[s]
            max_pro = curr_pro[s]

    return max_path


obs = ['normal', 'cold', 'dizzy']
result = Viterbit(obs, states, start_probability,
                  transition_probability, emission_probability)

print("Resultado:", result)
