#include("../src/DalitzPlot.jl")
#using .DalitzPlot, .DalitzPlot.Xs, .DalitzPlot.plot
using DalitzPlot, DalitzPlot.Xs, DalitzPlot.plot, DalitzPlot.QFT
#using DalitzPlot
using Test, ProgressBars
@testset "DalitzPlot Tests" begin

    function amp(tecm, kf, ch, para)

        # get kf as momenta in the center-of-mass ,
        #k1,k2,k3=getkf(kf)       
        #get kf as momenta in laboratory frame
        k1, k2, k3 = Xs.getkf(para.p, kf, ch)

        # Incoming particle momentum
        # Center-of-mass frame: p1 = [p 0.0 0.0 E1]
        #p1, p2 = pcm(tecm, ch.mi)
        # Laboratory frame
        p1, p2 = Xs.plab(para.p, ch.mi)

        #flux
        #flux factor for cross section
        fac = 1e9 / (4 * para.p * ch.mi[2] * (2 * pi)^5)

        k12 = k1 + k2
        s12 = k12*k12
        m = 3.2
        A = 1 / (s12 - m^2 + im * m * 0.1)

        total = abs2(A) * fac * 0.389379e-3

        return total
    end

    function main()
        ch = (mi=[1., 1.], mf=[1., 2., 3.],
            namei=["p^i_{1}", "p^i_{2}"], namef=["p^f_{1}", "p^f_{2}", "p^f_{3}"],
            amp=amp)

        p = 20.0

        function progress_callback(pb)
            ProgressBars.update(pb)  # 更新进度条
        end
        nevtot=Int64(1e7)
        pb = ProgressBar(1:nevtot)  # 创建进度条，范围从1到n
        callback = i -> progress_callback(pb)  # 创建回调函数，传入进度条对象
        res = Xs.Xsection(Xs.plab2pcm(p, ch.mi), ch, callback, axes=[23, 21], nevtot=nevtot, Nbin=1000, para=(p=p, l=1.0))
        @show Xs.plab2pcm(p, ch.mi), res.cs0
        plot.plotD(res)

    end

    main()

end





