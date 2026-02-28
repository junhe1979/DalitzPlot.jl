# using DalitzPlot, DalitzPlot.Xs, DalitzPlot.PLOT, DalitzPlot.FR
using DalitzPlot, DalitzPlot.Xs, DalitzPlot.PLOT, DalitzPlot.FR
using Test, ProgressBars

@testset "DalitzPlot Tests" begin

    # Amplitude function for the process
    function amps(tecm, kf, ch, para, p0)
        # Get kf as momenta in the laboratory frame
        k1, k2, k3 = Xs.getkf(para.p, kf, ch)

        # Incoming particle momentum in laboratory frame
        p1, p2 = Xs.plab(para.p, ch.mi)

        # Flux factor for cross section calculation
        fac = 1e9 / (4 * para.p * ch.mi[2] * (2 * pi)^5)

        k12 = k1 + k2
        s12 = k12 * k12
        m = 6.0
        A = 1 / (s12 - m^2 + im * m * 0.1)

        total = abs2(A) * fac * 0.389379e-3

        return total
    end

    function main()
        Ecm = 20.0
        nevtot = Int64(1e7)

        # Progress bar update callback
        function progress_callback(pb)
            ProgressBars.update(pb)  # Update progress bar
        end

        # Dummy amplitude function
        amps1(tecm, kf, proc, para, p0) = 1.

        # Two-particle channel definition
        proc = (pf=["p1", "p2"],
            mi=[1.0, 1.0], mf=[0.0 for i in 1:2],
            namei=["p^i_{1}", "p^i_{2}"], namef=["p_{1}", "p_{2}"],
            amps=amps1)

        pb = ProgressBar(1:nevtot)  # Create progress bar, range from 1 to nevtot
        callback = i -> progress_callback(pb)  # Create callback function, passing progress bar object
        res = Xs.Xsection(Ecm, proc, callback, axes=["p1:p2"], nevtot=nevtot, Nbin=1000,
            para=(p=Ecm, l=1.0), stype=2)
        @show Ecm, res.cs0 / (pi / 2.)

        # Three-particle channel definition (final state masses 0)
        proc = (pf=["p1", "p2", "p3"],
            mi=[1.0, 1.0], mf=[0.0 for i in 1:3],
            namei=["p^i_{1}", "p^i_{2}"], namef=["p_{1}", "p_{2}", "p_{3}"],
            amps=amps1)

        pb = ProgressBar(1:nevtot)  # Create progress bar, range from 1 to nevtot
        callback = i -> progress_callback(pb)  # Create callback function, passing progress bar object
        res = Xs.Xsection(Ecm, proc, callback, axes=["p3:p2", "p1:p2"], nevtot=nevtot, Nbin=1000,
            para=(p=Ecm, l=1.0), stype=2)

         @show Ecm, res.cs0 / (pi^2*Ecm^2/ 8.)

        # Three-particle channel definition (final state masses 2.0)
        proc = (pf=[:"p1", "p2", "p3"],
            mi=[1.0, 1.0], mf=[2.0 for i in 1:3],
            namei=["p^i_{1}", "p^i_{2}"], namef=["p_{1}", "p_{2}", "p_{3}"],
            amps=amps)
        pb = ProgressBar(1:nevtot)  # Create progress bar, range from 1 to nevtot
        callback = i -> progress_callback(pb)  # Create callback function, passing progress bar object
        res = Xs.Xsection(10.0, proc, callback, axes=["p2:p3", "p1:p2", "p1:p3"], nevtot=nevtot, Nbin=1000,
            para=(p=8000.0, l=1.0), stype=2)
        PLOT.plotD(res,filename="DP.png")
    end

    main()
end
