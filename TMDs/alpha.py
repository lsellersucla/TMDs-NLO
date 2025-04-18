import apfel

class ALPHAS:
    def __init__(self, order):
#        apfel.EnableWelcomeMessage(False)
        apfel.SetPerturbativeOrder(order)
        apfel.InitializeAPFEL()

    def alfs(self,Q):
        if Q < 0.5:
            Q = 0.5
        return apfel.AlphaQCD(Q)



