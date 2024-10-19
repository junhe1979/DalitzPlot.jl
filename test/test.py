#!/usr/bin/env python
# coding=utf-8
from juliacall import Main as jl
jl.seval("using StaticArrays")
jl.seval("using DalitzPlot")
jl.seval("using DalitzPlot.FR")
jl.seval("using DalitzPlot.GEN")
jl.seval("using DalitzPlot.Xs")
jl.seval("using DalitzPlot.qBSE")
jl.seval("using DalitzPlot.plot")

def amp(tecm, kf, ch, para):
    return 1.

# Main function to run the process
def main():
    # Define channel parameters
    ch = jl. NamedTuple{(:mi, :mf, :namei, :namef, :amp)}(mi=[1.0, 1.0],
        mf=[1.0, 2.0, 3.0],
        namei=["p^i_{1}", "p^i_{2}"],
        namef=["p^f_{1}", "p^f_{2}", "p^f_{3}"],
        amp=amp  # Assume `amp` is defined elsewhere
    )

    # Define other parameters
    p = 20.0
    nevtot = int(1e7)

   
    callback = lambda i: i

    # Call Xsection method (assuming a similar function in Python)
    res = jl.Xs.Xsection(jl.Xs.plab2pcm(p, ch['mi']), ch, callback, axes=[23, 21], nevtot=nevtot, Nbin=1000, para={'p': p, 'l': 1.0})

    # Show cross section results
    print(f"plab2pcm: {jl.Xs.plab2pcm(p, ch['mi'])}, cs0: {res['cs0']}")

    # Plot the results (assuming plot.plotD works similarly)
    jl.plot.plotD(res)

# Run the main function
if __name__ == '__main__':
    main()