#*******************************************************************************************
# Code for  Qusipotential Bethe-Salpter equation
# by Jun He @NNU   July 2024  
# Version 1: 29 July 2024 
# Version 2: 1 Sept 2024
#*******************************************************************************************
#using packages
using Distributed, ProgressBars #多进程并行包，进度条包
addprocs(5; exeflags="--project") #设定进程数
@everywhere using LinearAlgebra, StaticArrays, Printf #线性代数包，静态数组包，格式化输出包
#@everywhere include("../src/DalitzPlot.jl") #引入DalitzPlot包的路径
@everywhere using DalitzPlot, DalitzPlot.qBSE, DalitzPlot.FR
@everywhere qBSE.particles!(qBSE.p, qBSE.pkey,"data/particles.txt")  #引入粒子数据如质量量子数
#*******************************************************************************************
# 相互作用势。其中的例子为计算L1405和Zc3900。
@everywhere function fV(k, l, CB, IA, ichi, ichf) #
    if IA[ichi, ichf].Project == "L1405"
        U2i = FR.U(k.i2, l.i2)  # 求初态旋量
        U2f = FR.U(k.f2, l.f2, bar=true)  # # 求末态旋量

        GSk1 = FR.GS(k.i1 + k.f1)
        LLi, LLf = CB[ichi].cutoff, CB[ichf].cutoff #cutoff
        FF = qBSE.fPropFF(k, 70, (LLi + LLf) / 2.0, LLi, LLf, lregu=1, lFFex=0) #  求传播子和形状因子
        fV = -1.0 / 4.0 / (0.85 * 0.093)^2 * FF * (U2f * GSk1 * U2i) * IA[ichi, ichf].CC[1] #势
 

    elseif IA[ichi, ichf].Project == "Zc3900"
        beta, gV, fpi, gsigma, gcc, lambda, g2 = 0.90, 5.8, 0.132, 0.7614, 0.59, 0.56, 1.1618 #物理参数
        mi1, mi2 = CB[ichi].m[1], CB[ichi].m[2]  #质量
        mf1, mf2 = CB[ichf].m[1], CB[ichf].m[2]
        if CB[ichf].J[1] == 1 #如果粒子为矢量粒子，则给出极化矢量
            epsf1 = FR.eps(k.f1, l.f1, star=true)
        else
            epsf1 = SVector{5,ComplexF64}(0.0, 0.0, 0.0, 0.0, 0.0) #如果不是，这极化矢量设为0
        end
        if CB[ichi].J[1] == 1
            epsi1 = FR.eps(k.i1, l.i1)
        else
            epsi1 = SVector{5,ComplexF64}(0.0, 0.0, 0.0, 0.0, 0.0)
        end
        if CB[ichf].J[2] == 1
            epsf2 = FR.eps(k.f2, l.f2, star=true)
        else
            epsf2 = SVector{5,ComplexF64}(0.0, 0.0, 0.0, 0.0, 0.0)
        end
        if CB[ichi].J[2] == 1
            epsi2 = FR.eps(k.i2, l.i2)
        else
            epsi2 = SVector{5,ComplexF64}(0.0, 0.0, 0.0, 0.0, 0.0)
        end

        KerV = 0.0
        for le in 1:IA[ichi, ichf].Nex #遍历所有交换粒子
            name_ex = IA[ichi, ichf].name_ex[le] #交换粒子的名字
            key_ex = qBSE.pkey[name_ex] #交换粒子的指标
            J_ex, Jh_ex, m_ex = qBSE.p[key_ex].J, qBSE.p[key_ex].Jh, qBSE.p[key_ex].m #交换粒子的量子数及质量
            if IA[ichi, ichf].dc[le] == 1 #direct图
                name_ex = name_ex * "" #名字直接用交换粒子名字
                name1 = CB[ichi].name[1] * '-' * CB[ichf].name[1] #定义顶点1的名字
                name2 = CB[ichi].name[2] * '-' * CB[ichf].name[2] #定义顶点2的名字
                k.q = k.f2 - k.i2 #定义交换动量
                k.q2 = -abs(k.q * k.q)
            elseif IA[ichi, ichf].dc[le] == 2 #cross图
                name_ex = name_ex * "&c" #加&c
                name1 = CB[ichi].name[1] * '-' * CB[ichf].name[2] #定义顶点1的名字
                name2 = CB[ichi].name[2] * '-' * CB[ichf].name[1] #定义顶点2的名字
                k.q = k.f1 - k.i2
                k.q2 = -abs(k.q * k.q)
            end
            #以下是给出具体的顶点的费曼规则
            if name1 == "D-D" || name1 == "B-B"
                if name_ex in ["rho", "omega", "phi"]
                    Ver1 = im * beta * gV * sqrt(2.0) * (k.i1 + k.f1) / 2.0
                elseif name_ex in ["rho&c", "omega&c", "phi&c"]
                    Ver1 = im * beta * gV * sqrt(2.0) * (k.i1 + k.f2) / 2.0
                elseif name_ex == "Jpsi"
                    Ver1 = 2.0im * g2 * sqrt(m_ex) * mi1 * (k.i1 + k.f1)
                elseif name_ex == "Jpsi&c"
                    Ver1 = 2.0im * g2 * sqrt(m_ex) * mi1 * (k.i1 + k.f2)
                elseif name_ex == "Upsilon"
                    Ver1 = 2.0im * 2.0 * 0.40729 * sqrt(m_ex) * mi1 * (k.i1 + k.f1)
                elseif name_ex == "Upsilon&c"
                    Ver1 = 2.0im * 0.40729 * sqrt(m_ex) * mi1 * (k.i1 + k.f2)
                elseif name_ex == "sig"
                    Ver1 = -2.0im * gsigma * mi1
                elseif name_ex == "sig&c"
                    Ver1 = -2.0im * gsigma * mi1
                else
                    println("no such exchanged meson")
                end
            end
            if name2 == "D-D" || name2 == "B-B"
                if name_ex in ["rho", "omega", "phi"]
                    Ver2 = -1.0im * beta * gV * sqrt(2.0) * (k.i2 + k.f2) / 2.0
                elseif name_ex in ["rho&c", "omega&c", "phi&c"]
                    Ver2 = -1.0im * beta * gV * sqrt(2.0) * (k.i2 + k.f1) / 2.0
                elseif name_ex == "Jpsi"
                    Ver2 = -2.0im * g2 * sqrt(m_ex) * mi2 * (k.i2 + k.f2)
                elseif name_ex == "Jpsi&c"
                    Ver2 = -2.0im * g2 * sqrt(m_ex) * mi2 * (k.i2 + k.f1)
                elseif name_ex == "Upsilon"
                    Ver2 = -2.0im * 0.40729 * sqrt(m_ex) * mi2 * (k.i2 + k.f2)
                elseif name_ex == "Upsilon&c"
                    Ver2 = -2.0im * 0.40729 * sqrt(m_ex) * mi2 * (k.i2 + k.f1)
                elseif name_ex == "sig"
                    Ver2 = -2.0im * gsigma * mi2
                elseif name_ex == "sig&c"
                    Ver2 = -2.0im * gsigma * mi2
                else
                    println("no such exchanged meson")
                end
            end
            if name1 == "DA-DA" || name1 == "BA-BA"
                if name_ex in ["rho", "omega", "phi"]
                    Ver1 = -1.0im * beta * gV * sqrt(2.0) * (epsf1 * epsi1) * (k.i1 + k.f1) / 2.0 -
                           2.0im * sqrt(2.0) * lambda * gV * mi1 * (-(epsf1 * k.q) * epsi1 + (epsi1 * k.q) * epsf1)
                elseif name_ex in ["rho&c", "omega&c", "phi&c"]
                    Ver1 = -1.0im * beta * gV * sqrt(2.0) * (epsf2 * epsi1) * (k.i1 + k.f2) / 2.0 -
                           2.0im * sqrt(2.0) * lambda * gV * mi1 * (-(epsf2 * k.q) * epsi1 + (epsi1 * k.q) * epsf2)
                elseif name_ex == "Jpsi"
                    Ver1 = -2.0im * g2 * sqrt(m_ex) * mi1 * (-(epsi1 * (k.i1 + k.f1)) * epsf1 - (epsf1 * (k.i1 + k.f1)) * epsi1 + (epsf1 * epsi1) * (k.i1 + k.f1))
                elseif name_ex == "Jpsi&c"
                    Ver1 = -2.0im * g2 * sqrt(m_ex) * mi1 * (-(epsi1 * (k.i1 + k.f2)) * epsf2 - (epsf2 * (k.i1 + k.f2)) * epsi1 + (epsf2 * epsi1) * (k.i1 + k.f2))
                elseif name_ex == "Upsilon"
                    Ver1 = -2.0im * 0.40729 * sqrt(m_ex) * mi1 * (-(epsi1 * (k.i1 + k.f1)) * epsf1 - (epsf1 * (k.i1 + k.f1)) * epsi1 + (epsf1 * epsi1) * (k.i1 + k.f1))
                elseif name_ex == "Upsilon&c"
                    Ver1 = -2.0im * 0.40729 * sqrt(m_ex) * mi1 * (-(epsi1 * (k.i1 + k.f2)) * epsf2 - (epsf2 * (k.i1 + k.f2)) * epsi1 + (epsf2 * epsi1) * (k.i1 + k.f2))
                elseif name_ex in ["pi", "eta", "etap"]
                    Ver1 = -2.0im * gcc / fpi * FR.LC((k.f1 + k.i1) / 2.0, epsi1, -k.q, epsf1)
                elseif name_ex in ["pi&c", "eta&c", "etap&c"]
                    Ver1 = -2.0im * gcc / fpi * FR.LC((k.f2 + k.i1) / 2.0, epsi1, -k.q, epsf2)
                elseif name_ex == "sig"
                    Ver1 = 2.0im * gsigma * mi1 * (epsf1 * epsi1)
                elseif name_ex == "sig&c"
                    Ver1 = 2.0im * gsigma * mi1 * (epsf2 * epsi1)
                else
                    println("no such exchanged meson")
                end
            end
            if name2 == "DA-DA" || name2 == "BA-BA"
                if name_ex in ["rho", "omega", "phi"]
                    Ver2 = 1.0im * sqrt(2.0) * beta * gV * (epsf2 * epsi2) * (k.i2 + k.f2) / 2.0 -
                           2.0im * sqrt(2.0) * lambda * gV * mi2 * ((epsf2 * k.q) * epsi2 - (epsi2 * k.q) * epsf2)
                elseif name_ex in ["rho&c", "omega&c", "phi&c"]
                    Ver2 = 1.0im * sqrt(2.0) * beta * gV * (epsf1 * epsi2) * (k.i2 + k.f1) / 2.0 -
                           2.0im * sqrt(2.0) * lambda * gV * mi2 * ((epsf1 * k.q) * epsi2 - (epsi2 * k.q) * epsf1)
                elseif name_ex == "Jpsi"
                    Ver2 = 2.0im * g2 * sqrt(m_ex) * mi2 * (-(epsi2 * (k.i2 + k.f2)) * epsf2 - (epsf2 * (k.i2 + k.f2)) * epsi2 + (epsf2 * epsi2) * (k.i2 + k.f2))
                elseif name_ex == "Jpsi&c"
                    Ver2 = 2.0im * g2 * sqrt(m_ex) * mi2 * (-(epsi2 * (k.i2 + k.f1)) * epsf1 - (epsf1 * (k.i2 + k.f1)) * epsi2 + (epsf1 * epsi2) * (k.i2 + k.f1))
                elseif name_ex == "Upsilon"
                    Ver2 = 2.0im * 0.40729 * sqrt(m_ex) * mi2 * (-(epsi2 * (k.i2 + k.f2)) * epsf2 - (epsf2 * (k.i2 + k.f2)) * epsi2 + (epsf2 * epsi2) * (k.i2 + k.f2))
                elseif name_ex == "Upsilon&c"
                    Ver2 = 2.0im * 0.40729 * sqrt(m_ex) * mi2 * (-(epsi2 * (k.i2 + k.f1)) * epsf1 - (epsf1 * (k.i2 + k.f1)) * epsi2 + (epsf1 * epsi2) * (k.i2 + k.f1))
                elseif name_ex in ["pi", "eta", "etap"]
                    Ver2 = -2.0im * gcc / fpi * FR.LC((k.f2 + k.i2) / 2.0, epsi2, k.q, epsf2)
                elseif name_ex in ["pi&c", "eta&c", "etap&c"]
                    Ver2 = -2.0im * gcc / fpi * FR.LC((k.f1 + k.i2) / 2.0, epsi2, k.q, epsf1)
                elseif name_ex == "sig"
                    Ver2 = 2.0im * gsigma * mi2 * (epsf2 * epsi2)
                elseif name_ex == "sig&c"
                    Ver2 = 2.0im * gsigma * mi2 * (epsf1 * epsi2)
                else
                    println("no such exchanged meson")
                end
            end
            if name1 == "DA-D" || name1 == "BA-B"
                if name_ex in ["rho", "omega", "phi"]
                    tempV2 = -FR.LC((k.i1 + k.f1) / 2.0, -k.q, epsi1)
                    Ver1 = -2.0 * sqrt(2) * lambda * gV * tempV2
                elseif name_ex in ["rho&c", "omega&c", "phi&c"]
                    tempV2 = -FR.LC((k.i1 + k.f2) / 2.0, -k.q, epsi1)
                    Ver1 = -2.0 * sqrt(2) * lambda * gV * tempV2
                elseif name_ex == "Jpsi"
                    tempV2 = -FR.LC(k.i1 + k.f1, -k.q, epsi1)
                    tempV = sqrt(2) * lambda * gV * tempV2
                    Ver1 = 2.0 * g2 * sqrt(m_ex) * tempV2
                elseif name_ex == "Jpsi&c"
                    tempV2 = -FR.LC(k.i1 + k.f2, -k.q, epsi1)
                    tempV = sqrt(2) * lambda * gV * tempV2
                    Ver1 = 2.0 * g2 * sqrt(m_ex) * tempV2
                elseif name_ex == "Upsilon"
                    tempV2 = -FR.LC(k.i1 + k.f1, -k.q, epsi1)
                    tempV = sqrt(2) * lambda * gV * tempV2
                    Ver1 = 2.0 * 0.40729 * sqrt(m_ex) * tempV2
                elseif name_ex == "Upsilon&c"
                    tempV2 = -FR.LC(k.i1 + k.f2, -k.q, epsi1)
                    tempV = sqrt(2) * lambda * gV * tempV2
                    Ver1 = 2.0 * 0.40729 * sqrt(m_ex) * tempV2
                elseif name_ex in ["pi", "eta", "etap"]
                    Ver1 = -2.0 * gcc / fpi * sqrt(mi1 * mf1) * (epsi1 * k.q)
                elseif name_ex in ["pi&c", "eta&c", "etap&c"]
                    Ver1 = -2.0 * gcc / fpi * sqrt(mi1 * mf2) * (epsi1 * k.q)
                else
                    println("no such exchanged meson")
                end
            end
            if name2 == "D-DA" || name2 == "B-BA"
                if name_ex in ["rho", "omega", "phi"]
                    tempV2 = -FR.LC((k.i2 + k.f2) / 2.0, k.q, epsf2)
                    Ver2 = -2.0 * sqrt(2) * lambda * gV * tempV2
                elseif name_ex in ["rho&c", "omega&c", "phi&c"]
                    tempV2 = -FR.LC((k.i2 + k.f1) / 2.0, k.q, epsf1)
                    Ver2 = -2.0 * sqrt(2) * lambda * gV * tempV2
                elseif name_ex == "Jpsi"
                    tempV2 = -FR.LC(k.i2 + k.f2, k.q, epsf2)
                    Ver2 = 2.0 * g2 * sqrt(p[p_Jpsi].m) * tempV2
                elseif name_ex == "Jpsi&c"
                    tempV2 = -FR.LC(k.i2 + k.f1, k.q, epsf1)
                    Ver2 = 2.0 * g2 * sqrt(p[p_Jpsi].m) * tempV2
                elseif name_ex == "Upsilon"
                    tempV2 = -FR.LC(k.i2 + k.f2, k.q, epsf2)
                    Ver2 = 2.0 * 0.40729 * sqrt(m_ex) * tempV2
                elseif name_ex == "Upsilon&c"
                    tempV2 = -FR.LC(k.i2 + k.f1, k.q, epsf1)
                    Ver2 = 2.0 * 0.40729 * sqrt(m_ex) * tempV2
                elseif name_ex in ["pi", "eta", "etap"]
                    Ver2 = -2.0 * gcc / fpi * sqrt(mi2 * mf2) * (epsf2 * k.q)
                elseif name_ex in ["pi&c", "eta&c", "etap&c"]
                    Ver2 = -2.0 * gcc / fpi * sqrt(mi2 * mf1) * (epsf1 * k.q)
                else
                    println("no such exchanged meson")
                end
            end
            if name1 == "D-DA" || name1 == "B-BA"
                if name_ex in ["rho", "omega", "phi"]
                    tempV2 = -FR.LC((k.i1 + k.f1) / 2, -k.q, epsf1)
                    Ver1 = -2.0 * sqrt(2) * lambda * gV * tempV2
                elseif name_ex in ["rho&c", "omega&c", "phi&c"]
                    tempV2 = -FR.LC((k.i1 + k.f2) / 2, -k.q, epsf2)
                    Ver1 = -2.0 * sqrt(2) * lambda * gV * tempV2
                elseif name_ex == "Jpsi"
                    tempV2 = -FR.LC(k.i1 + k.f1, -k.q, epsf1)
                    Ver1 = 2.0 * g2 * sqrt(m_ex) * tempV2
                elseif name_ex == "Jpsi&c"
                    tempV2 = -FR.LC(k.i1 + k.f2, -k.q, epsf2)
                    Ver1 = 2.0 * g2 * sqrt(m_ex) * tempV2
                elseif name_ex == "Upsilon"
                    tempV2 = -FR.LC(k.i1 + k.f1, -k.q, epsf1)
                    Ver1 = 2.0 * 0.40729 * sqrt(m_ex) * tempV2
                elseif name_ex == "Upsilon&c"
                    tempV2 = -FR.LC(k.i1 + k.f2, -k.q, epsf2)
                    Ver1 = 2.0 * 0.40729 * sqrt(m_ex) * tempV2
                elseif name_ex in ["pi", "eta", "etap"]
                    Ver1 = -2.0 * gcc / fpi * sqrt(mi1 * mf1) * (epsf1 * k.q)
                elseif name_ex in ["pi&c", "eta&c", "etap&c"]
                    Ver1 = -2.0 * gcc / fpi * sqrt(mi1 * mf2) * (epsf2 * k.q)
                else
                    println("no such exchanged meson")
                end
            end
            if name2 == "DA-D" || name2 == "BA-B"
                if name_ex in ["rho", "omega", "phi"]
                    tempV2 = -FR.LC((k.i2 + k.f2) / 2.0, k.q, epsi2)
                    Ver2 = -2.0 * sqrt(2.0) * lambda * gV * tempV2
                elseif name_ex in ["rho&c", "omega&c", "phi&c"]
                    tempV2 = -FR.LC((k.i2 + k.f1) / 2, k.q, epsi2)
                    Ver2 = -2.0 * sqrt(2.0) * lambda * gV * tempV2
                elseif name_ex == "Jpsi"
                    tempV2 = -FR.LC(k.i2 + k.f2, k.q, epsi2)
                    Ver2 = 2.0 * g2 * sqrt(m_ex) * tempV2
                elseif name_ex == "Jpsi&c"
                    tempV2 = -FR.LC(k.i2 + k.f1, k.q, epsi2)
                    Ver2 = 2.0 * g2 * sqrt(m_ex) * tempV2
                elseif name_ex == "Upsilon"
                    tempV2 = -FR.LC(k.i2 + k.f2, k.q, epsi2)
                    Ver2 = 2.0 * 0.40729 * sqrt(m_ex) * tempV2
                elseif name_ex == "Upsilon&c"
                    tempV2 = -FR.LC(k.i2 + k.f1, k.q, epsi2)
                    Ver2 = 2.0 * 0.40729 * sqrt(m_ex) * tempV2
                elseif name_ex in ["pi", "eta", "etap"]
                    Ver2 = -2.0 * gcc / fpi * sqrt(mi2 * mf2) * (epsi2 * k.q)
                elseif name_ex in ["pi&c", "eta&c", "etap&c"]
                    Ver2 = -2.0 * gcc / fpi * sqrt(mi2 * mf1) * (epsi2 * k.q)
                else
                    println("no such exchanged meson")
                end
            end
            LLi, LLf = CB[ichi].cutoff, CB[ichf].cutoff
            FF = qBSE.fPropFF(k, key_ex, (LLi + LLf) / 2.0, LLi, LLf, lFFex=3)  
            if J_ex == 0 #标量赝标交换
                KerV -= Ver1 * Ver2 * IA[ichi, ichf].CC[le] * FF #由顶点求势
            end
            if J_ex == 1 && Jh_ex == 1 # 矢量交换
                KerV += ((Ver1 * Ver2) - (Ver1 * k.q) * (Ver2 * k.q) / m_ex^2) * IA[ichi, ichf].CC[le] * FF
            end
        end
        fV = KerV
    end
    return fV
end
#*******************************************************************************************
# 并行运算的内容，主要是计算需要的量，比如这里是得到log|1-VG|以寻找极点。
@everywhere function res(Range, iER, qn, CB, IH, IA, fV)
    resEct = ComplexF64[] #设置保存复的总能量Ec=ER+EI*im的数组
    reslogt = Float64[] #设置保存log|1-VG|的数组
    ER = Range.ERmax - iER * (Range.ERmax - Range.ERmin) / Range.NER #计算ER
    for iEI in -Range.NEI:Range.NEI #虚部部分循环
        EI = iEI * Range.EIt / Range.NEI #计算虚部
        Ec = ER + EI * im
        Vc, Gc, II = qBSE.srAB(Ec, qn, CB, IH, IA, fV,lRm=1) # Calculate the V, G, and unit matrix II
        VGI = II - Vc * Gc
        detVGI = det(VGI)   # Compute determinant of (1 - VGc)，调用LinearAlgebra包det函数
        logdetVGI = log(abs(detVGI)^2)
        push!(resEct, Ec)
        push!(reslogt, logdetVGI)
    end
    return resEct, reslogt
end
#*******************************************************************************************
function main() #主函数，主要的计算流程
    #----------------------------------------------------------------
    dashline = repeat('-', 100) 
    println(dashline)
    println("Start project: ")
    #############################################################################
    #输入系统信息
    Project = "L1405" #选择计算L1405还是Zc3900
    if Project == "L1405"  
        qn = (I=1, Ih=2, J=1, Jh=2, P=-1, C=-1) #体系量子数I,Ih,J,Jh,P,C 角动量J/Jh，同位旋为I/Ih            
        Range = (ERmin=1.2, ERmax=1.6, NER=200, EIt=0.200, NEI=20) #计算的范围及精度
        channels = ( #系统包含的道
            ("K", "N", qBSE.p[qBSE.pkey["K"]].m + 0.22 * 1.63 * 0.5), #道的粒子及cutoff
            ("pi", "Sigma", qBSE.p[qBSE.pkey["pi"]].m + 0.22 * 1.63 * 0.5))
        #被保存到IA， 其中标签转化为[chi,chf], 后面依次为Nex,ex,CC
        CC = Dict(  #flavor factors
            "K:N-->K:N" => ([["V", 1], 3.0],), #过程对应的交换，这里为contact图，以V表示，最后为cutoff
            "K:N-->pi:Sigma" => ([["V", 1], -sqrt(1.5)],),
            "pi:Sigma-->K:N" => ([["V", 1], -sqrt(1.5)],),
            "pi:Sigma-->pi:Sigma" => ([["V", 1], 4.0],))
    end
    #
    if Project == "Zc3900"
        qn = (I=1, Ih=1, J=1, Jh=1, P=1, C=-1) #I,Ih,J,Jh,P,C               
        Range = (ERmin=3.850, ERmax=3.880, NER=10, EIt=0.040, NEI=10)
        channels = (("D", "DA", 3.3),)
        CC = Dict(
            "D:DA-->D:DA" => (
                [["rho", 1], -0.5], #1,2为direct或cross图
                [["omega", 1], 0.5],
                [["sig", 1], 1.0],
                [["Jpsi", 1], 1.0],
                [["rho", 2], -0.5 * (-qn.C)],
                [["omega", 2], 0.5 * (-qn.C)],
                [["pi", 2], -0.5 * (-qn.C)],
                [["eta", 2], 1.0 / 6.0 * (-qn.C)],
                [["Jpsi", 2], 1.0]),
            "DA:DA-->DA:DA" => ( # 考虑DADA道是有用
                [["rho", 1], -0.5],
                [["omega", 1], 0.5],
                [["pi", 1], -0.5],
                [["eta", 1], 1.0 / 6.0],
                [["sig", 1], 1.0],
                [["Jpsi", 1], 1.0])
        )
    end
    #----------------------------------------------------------------
    CB, IH, IA = qBSE.Independent_amp(Project, channels, CC, qn)#生成包含道，独立振幅和相互作用信息的struct实例
    #----------------------------------------------------------------
    println(dashline)
    println(dashline)
    println("I(J,P)=$(qn.I)($(qn.J),$(qn.P))  ")
    for ih in eachindex(IH)
        println("")
        println("$(CB[IH[ih].ich].name[3]) cutoff = $(CB[IH[ih].ich].cutoff) GeV: $(IH[ih].hel[1]) / $(IH[ih].helh[1])  $(IH[ih].hel[2]) / $(IH[ih].helh[2])")
    end
    Nc = length(CB)
    println(dashline)
    println("channel        ", join([CB[i].name[3] for i in 1:Nc], " "))
    for i1 in 1:Nc
        println(CB[i1].name[3], join([IA[i1, i].Nex for i in 1:Nc], " "))
    end
    println(dashline)
    println("ER=$(Range.ERmax) to $(Range.ERmin) GeV, NER=$(Range.NER)  ;   EI=$(Range.EIt*1e3) MeV, NEI=$(Range.NEI)")
    println(dashline)
    labelname = [CB[itemp].name[1] for itemp in 1:Nc]
    println("E  D= :: ", join(labelname[1:Nc], " "))
    #############################################################################
    #计算
    if isfile("outputjl.txt")    #清除旧数据
        run(`rm outputjl.txt`)
    end
    resEc = ComplexF64[]
    reslog = Float64[]
    Nprocs = nprocs() - 1
    for iER0 in ProgressBar(0:Range.NER/Nprocs) #实部循环
        #results=[res(Range, iER0*Nprocs,qn, CB, IH, IA, fV)]
        results = pmap(i -> res(Range, iER0 * Nprocs + i - 1, qn, CB, IH, IA, fV), 1:Nprocs) #多进程
        for res0 in results
            append!(resEc, res0[1])
            append!(reslog, res0[2])
        end
    end
    #--------------------------------------------------------------------
    # 寻找极点
    Ampmin, Ampminx, Ampminy = 0.0, 0.0, 0.0
    for i in eachindex(resEc)
        Ec = resEc[i]
        logdetVGI = reslog[i]
        if abs(imag(Ec)) < 0.000001
            #@printf("Output: Ec=%.4f, logD=%.4f\n", real(Ec), logdetVGI)
        end
        open("outputjl.txt", "a") do file
            write(file, @sprintf("%.4f, %.4f, %.4f\n", real(Ec), imag(Ec) * 1e3, logdetVGI))
        end
        if Ampmin > logdetVGI
            Ampmin = logdetVGI
            Ampminx = real(Ec)
            Ampminy = imag(Ec) * 1e3
        end
    end
    #----------------------------------------------------------------
    println("I(J,P)=$(qn.I)($(qn.J),$(qn.P))  pole= $Ampmin at $(Ampminx * 1e3), $Ampminy")
    println(dashline)
    println("END A LOOP")
    println(dashline)
    println(dashline)
    println("END PROGRAM")
    println(dashline)
end
@time main() #运行主程序

include("plot.jl") #gnuplot画图
#run(`gnuplot plot/plotL1405jl.gp`) #matplotlib画图