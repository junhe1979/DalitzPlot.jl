using StaticArrays
using Test

using DalitzPlot,StaticArrays

@testset "DalitzPlot.jl" begin
    
    function amp(kf, ch, para)

        #入射流
        #flux factor for cross section
        fac = 1 / (4 * para.pW * ch.p_i[2] * (2 * pi)^5)
        #fermion, aevrage on initial particele 
        fac = fac * 2.0 * ch.p_i[2] * 2.0 * ch.p_f[2] / 2.0
        
        #产生的kf为质心系动量，质心系->实验室系
        k1,k2,k3=kcm2klab(para.pW, kf,ch)
        
        #质心系下入射粒子动量
        p1,p2=plab(para.pW, ch.p_i)

       
        k12=k1+k2
        s12 = cdot(k12, k12)
        m=3.
        A = 1/ (s12 - m ^2 + im * m * 0.1)
    
        total = abs2(A) * fac
    
        return total
    end
    function main()
    ch = (p_i=[1., 1.], p_f=[1., 1., 1.],
        Lp_i=["p^i_{1}", "p^i_{2}"], Lp_f=["p^f_{1}", "p^f_{2}", "p^f_{3}"],
        amp=amp) #初末态粒子质量

    p= 10.0    
    res=Xsection(plab2pcm(p,ch.p_i), ch,nevtot=Int64(1e7),para=(pW=p, l=1.0),ProgressBars=true)
    @show p, res.cs0
    plotD(res,ch,axes=[1, 3])
    end 

    main()

end
