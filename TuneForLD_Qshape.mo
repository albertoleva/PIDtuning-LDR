package TuneForLD_Qshape
  package Constants
    constant Integer nClasses = 5 "No. of process classes";
    constant Integer nRules = 10 "No. of tuning rules";
    constant Integer pmminRp = 30 "default min pm for proposed rule";
    constant Boolean enforce_pmminRp = true "enforce the above";
    constant Integer MsVAP18_10 = 20 "Ms for VAP18 times 10 (e.g., 16 for 1.6)";
    constant Real NVAP18 = 20;
  end Constants;

  package Functions
    package Internal
      function factorial
        input Integer n;
        output Integer f;
      algorithm
        f := if n >= 1 then n * factorial(n - 1) else 1;
      end factorial;

      function beta1 "service function for analytical MoA"
        input Real p;
        input Real k;
        output Real y;
      algorithm
        y := exp(-(p + 1) * (p ^ 2 + 1) / p ^ k) - 1;
      end beta1;

      function PID_classical2isa "classical PID K(1+1/sTi)(1+sTd)/(1+sTd/N) to ISA"
        input Real Kc, Tic, Tdc, Nc "classical PID parameters";
        output Real Ki, Tii, Tdi, Ni "ISA PID parameters";
      algorithm
        Ki := ((Kc * Nc - Kc) * Tdc + Kc * Nc * Tic) / (Nc * Tic);
        Tii := ((Nc - 1) * Tdc + Nc * Tic) / Nc;
        Tdi := ((1 - Nc) * Tdc ^ 2 + (Nc ^ 2 - Nc) * Tdc * Tic) / (Nc ^ 2 * Tic + (Nc ^ 2 - Nc) * Tdc);
        Ni := ((1 - Nc) * Tdc + (Nc ^ 2 - Nc) * Tic) / (Nc * Tic + (Nc - 1) * Tdc);
      end PID_classical2isa;

      function vec2mat_byrow "fill i-th row of mat (ncols columns) with i-th element of vec"
        input Real[:] vec;
        input Integer ncols;
        output Real[size(vec, 1), ncols] mat;
      algorithm
        for i in 1:size(vec, 1) loop
          mat[i, :] := vec[i] * ones(ncols);
        end for;
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
      end vec2mat_byrow;

      function write_matrix_to_csv "write matrix data to fileName in csv form"
        input Real[:, :] data;
        input String fileName;
      protected
        String curLine;
        Integer nrows, ncols;
      algorithm
        Modelica.Utilities.Files.remove(fileName);
        nrows := size(data, 1);
        ncols := size(data, 2);
        for i in 1:nrows loop
          curLine := "";
          for j in 1:ncols - 1 loop
            curLine := curLine + String(data[i, j]) + ",";
          end for;
          curLine := curLine + String(data[i, ncols]);
          Modelica.Utilities.Streams.print(curLine, fileName);
        end for;
        Modelica.Utilities.Streams.close(fileName);
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
      end write_matrix_to_csv;
    end Internal;

    package AnalyticalMoA "analytical method of areas for process classes 1-5 in Astrom & Hagglund2000 plus nominal case (1,1 Pade FOPDT)"
    
      function MoA_Pnom "1,1 Pade FOPDT with norm delay p"
        input Real p;
        output Real mu, T, D;
      algorithm
        mu := 1;
        T := 1;
        D := p/(1-p);
      end MoA_Pnom;
    
      function MoA_P1 "FOPDT approx for 1/(1+s)^n"
        input Real p;
        output Real mu, T, D;
      protected
        Integer n;
      algorithm
        n := integer(p);
        mu := 1;
        T := exp(1 - n) * n ^ n / Functions.Internal.factorial(n - 1);
        D := n - T;
      end MoA_P1;

      function MoA_P2 "FOPDT approx for 1/((1+s)(1+ps)(1+p^2s)(1+p^3s)"
        input Real p;
        output Real mu, T, D;
      protected
        Real A0, A1, beta2;
      algorithm
        mu := 1;
        if p < 1 then
          A0 := (p + 1) * (p ^ 2 + 1);
          beta2 := (p - 1) ^ 3 * (p + 1);
          A1 := A0 - (Functions.Internal.beta1(p, 2) * p ^ 5 - Functions.Internal.beta1(p, 1) * p ^ 2) / beta2 + (Functions.Internal.beta1(p, 3) * p ^ 9 - Functions.Internal.beta1(p, 0)) / (beta2 * (p ^ 2 + p + 1));
        else
          A0 := 4;
          A1 := 128 / 3 * exp(-4);
        end if;
        T := exp(1) * A1 / mu;
        D := A0 / mu - T;
      end MoA_P2;

      function MoA_P3 "FOPDT approx for (1-p*s)/(1+s)^3"
        input Real p;
        output Real mu, T, D;
      algorithm
        mu := 1;
        T := 1 / 2 * exp(-(p + 2)) * (p + 3) ^ 3;
        D := p + 3 - 0.5 * exp(-(p + 2)) * (p + 3) ^ 3;
      end MoA_P3;

      function MoA_P4 "FOPDT approx for exp(-s)/(1+p*s)"
        input Real p;
        output Real mu, T, D;
      algorithm
        mu := 1;
        T := p;
        D := 1;
      end MoA_P4;

      function MoA_P5 "FOPDT approx for exp(-s)/(1+p*s)^2"
        input Real p;
        output Real mu, T, D;
      algorithm
        mu := 1;
        T := 4 * exp(-1) * p;
        D := 1 + 2 * p * (1 - 2 * exp(-1));
      end MoA_P5;

      function MoA_P6 "FOPDT approx for (p^2)/(1+s)/(s^2+2*xi*s+p^2), xi=0.5"
        input Real p;
        output Real mu, T, D;
      algorithm
        mu := 1;
        T := exp((-1 / 2) - 3 / (2 * p ^ 2)) * (p ^ 2 * sqrt(4 * p ^ 2 - 1) * exp(1 / 2 + 1 / (2 * p ^ 2)) + sqrt(4 * p ^ 2 - 1) * exp(1 + 1 / p ^ 2) * cos((1 + p ^ 2) * sqrt(4 * p ^ 2 - 1) / (2 * p ^ 2)) + exp(1 + 1 / p ^ 2) * sin((1 + p ^ 2) * sqrt(4 * p ^ 2 - 1) / (2 * p ^ 2))) / (p ^ 2 * sqrt(4 * p ^ 2 - 1));
        D := -exp((-1 / 2) - 3 / (2 * p ^ 2)) * (p ^ 2 * sqrt(4 * p ^ 2 - 1) * exp(1 / 2 + 1 / (2 * p ^ 2)) - sqrt(4 * p ^ 2 - 1) * exp(1 / 2 + 3 / (2 * p ^ 2)) - p ^ 2 * sqrt(4 * p ^ 2 - 1) * exp(1 / 2 + 3 / (2 * p ^ 2)) + sqrt(4 * p ^ 2 - 1) * exp(1 + 1 / p ^ 2) * cos((1 + p ^ 2) * sqrt(4 * p ^ 2 - 1) / (2 * p ^ 2)) + exp(1 + 1 / p ^ 2) * sin((1 + p ^ 2) * sqrt(4 * p ^ 2 - 1) / (2 * p ^ 2))) / (p ^ 2 * sqrt(4 * p ^ 2 - 1));
      end MoA_P6;

      function MoA_P_all "FOPDT approx for process of class c with parameter value p"
        input Integer c(min = 1, max = Constants.nClasses) "process class (1-Constants.nClasses)";
        input Real p "process param value";
        output Real mu, T, D;
      algorithm
        if c == 0 then
          (mu, T, D) := Functions.AnalyticalMoA.MoA_Pnom(p);
        elseif c == 1 then
          (mu, T, D) := Functions.AnalyticalMoA.MoA_P1(p);
        elseif c == 2 then
          (mu, T, D) := Functions.AnalyticalMoA.MoA_P2(p);
        elseif c == 3 then
          (mu, T, D) := Functions.AnalyticalMoA.MoA_P3(p);
        elseif c == 4 then
          (mu, T, D) := Functions.AnalyticalMoA.MoA_P4(p);
        elseif c == 5 then
          (mu, T, D) := Functions.AnalyticalMoA.MoA_P5(p);
        else
          (mu, T, D) := Functions.AnalyticalMoA.MoA_P6(p);
        end if;
      end MoA_P_all;
    end AnalyticalMoA;

    package TuningRules "tuning rules"
      function PID_KS88IAEr "Kaya & Scheib 1988, min IAE, regulatory, classical PID"
        input Real mu, T, D;
        output Real K, Ti, Td, N;
      algorithm
        K := 0.98 / mu * (T / D) ^ 0.76;
        Ti := T / 0.91 * (T / D) ^ 1.05;
        Td := 0.6 * T * (D / T) ^ 0.9;
        N := 10;
// as per authors' proposal
      end PID_KS88IAEr;

      function PID_TS95dir "Tsang & Rad 1995, direct synthesis, classical PID"
        input Real mu, T, D;
        output Real K, Ti, Td, N;
      algorithm
        K := 0.81 * T / (mu * D);
        Ti := T;
        Td := 0.5 * D;
        N := 5;
// as per authors' proposal
      end PID_TS95dir;

      function PID_LC04std "Leva & Colombo 2004, IMC std lambda, ISA PID"
        input Real mu, T, D;
        output Real K, Ti, Td, N;
      protected
        Real lambda;
      algorithm
        lambda := max(0.25 * D, 0.20 * T);
        Ti := T + D ^ 2 / (2 * (lambda + D));
        K := Ti / (mu * (lambda + D));
        N := T * (lambda + D) / (lambda * Ti) - 1;
        Td := lambda * D * N / (2 * (lambda + D));
      end PID_LC04std;

      function PID_LC04small "Leva & Colombo 2004, IMC small lambda, ISA PID"
        input Real mu, T, D;
        output Real K, Ti, Td, N;
      protected
        Real lambda;
      algorithm
        lambda := 0.1 * max(0.25 * D, 0.20 * T);
        Ti := T + D ^ 2 / (2 * (lambda + D));
        K := Ti / (mu * (lambda + D));
        N := T * (lambda + D) / (lambda * Ti) - 1;
        Td := lambda * D * N / (2 * (lambda + D));
      end PID_LC04small;

      function PID_AR10r "Arrieta & al. 2010, regulatory, ISA PID (N=20 as per ref)"
        input Real mu, T, D;
        output Real K, Ti, Td, N;
      protected
        Real tau, a1, a2, a3, b1, b2, b3;
      algorithm
        tau := D / T;
        if tau <= 1.00 then
          a1 := 1.473 "values for tau in suggested range 0.1-1";
          b1 := -0.970;
          a2 := 1.115;
          b2 := -0.753;
          a3 := 0.550;
          b3 := 0.948;
        else
          a1 := 1.524 "values for tau in suggested range 1.1-2";
          b1 := -0.735;
          a2 := 1.130;
          b2 := -0.641;
          a3 := 0.552;
          b3 := 0.851;
        end if;
        K := a1 / mu * tau ^ b1;
        Ti := T / (a2 * tau ^ b2);
        Td := T * a3 * tau ^ b3;
        N := 20 "as per authors' proposal";
      end PID_AR10r;

      function PID_AMIGO "AMIGO, ISA PID, N as per AH 2004"
        input Real mu, T, D;
        output Real K, Ti, Td, N;
      protected
        Real theta;
      algorithm
        theta := D / (T + D);
        K := 1 / mu * (0.2 + 0.45 * T / D);
        Ti := (0.4 * D + 0.8 * T) / (D + 0.1 * T) * D;
        Td := 0.5 * D * T / (0.3 * D + T);
        N := -((-180) - 265 * theta + 1025 * theta ^ 2 - 665 * theta ^ 3 + 85 * theta ^ 4) / ((-31.5 * theta ^ 4) + 161.2 * theta ^ 3 - 266.1 * theta ^ 2 + 130.4 * theta + 18.0);
// figure 11 in AH JPC 2004
      end PID_AMIGO;

      function PI_SIMC "SIMC, tight control (FOPDT->PI)"
        input Real mu, T, D;
        output Real K, Ti, Td, N;
      protected
        Real tauc;
      algorithm
        tauc := D;
        K := T / mu / (tauc + D);
        Ti := min(T, 4 * (tauc + D));
        Td := 0;
        N := 1;
      end PI_SIMC;

      function PID_VAP18NVAP "Vilanova, Arrieta & Ponsa, ISA trans. 2018"
        input Real mu, T, D;
        input Integer Ms_10 = Constants.MsVAP18_10;
        input Real NVAP = Constants.NVAP18;
        output Real K, Ti, Td, N;
      protected
        Real p1, p2, p3, p4, tauo, tauc, kp, taui, taud;
      algorithm
        if Ms_10 == 20 then
          p1 := 0.73;
          p2 := 0.26;
          p3 := -0.72;
          p4 := -1.24;
        elseif Ms_10 == 18 then
          p1 := 0.75;
          p2 := 0.28;
          p3 := -0.73;
          p4 := -1.30;
        elseif Ms_10 == 16 then
          p1 := 0.76;
          p2 := 0.30;
          p3 := -0.73;
          p4 := -1.40;
        else
          p1 := 0;
          p2 := 0;
          p3 := 0;
          p4 := 0;
        end if;
        tauo := D / T;
        tauc := p1 * exp(p2 * tauo) + p3 * exp(p4 * tauo);
        kp := (tauo * (2 + tauo / 2) * (3 * tauc + tauo / 2) - 2 * tauc ^ 3 - 3 * tauc ^ 2 * tauo) / 2 / (tauc + tauo / 2) ^ 3;
        taui := (tauo * (2 + tauo / 2) * (3 * tauc + tauo / 2) - 2 * tauc ^ 3 - 3 * tauc ^ 2 * tauo) / 2 / tauo / (1 + tauo);
        taud := (3 * tauc ^ 2 * tauo + tauo ^ 2 / 2 * (3 * tauc + tauo / 2) - 2 * (1 + tauo) * tauc ^ 3) / (tauo * (2 + tauo / 2) * (3 * tauc + tauo / 2) - 2 * tauc ^ 3 - 3 * tauc ^ 2 * tauo);
        K := kp / mu;
        Ti := taui * T;
        Td := taud * T;
        N := NVAP;
      end PID_VAP18NVAP;

      function PID_VAP18NpropCinf "Vilanova, Arrieta & Ponsa, ISA trans. 2018"
        input Real mu, T, D;
        input Integer Ms_10 = Constants.MsVAP18_10;
        output Real K, Ti, Td, N;
      protected
        Real p1, p2, p3, p4, tauo, tauc, kp, taui, taud;
        Real kk, tti, ttd, nn;
      algorithm
        (kk, tti, ttd, nn) := Functions.TuningRules.PID_prop(mu, T, D);
        if Ms_10 == 20 then
          p1 := 0.73;
          p2 := 0.26;
          p3 := -0.72;
          p4 := -1.24;
        elseif Ms_10 == 18 then
          p1 := 0.75;
          p2 := 0.28;
          p3 := -0.73;
          p4 := -1.30;
        elseif Ms_10 == 16 then
          p1 := 0.76;
          p2 := 0.30;
          p3 := -0.73;
          p4 := -1.40;
        else
          p1 := 0;
          p2 := 0;
          p3 := 0;
          p4 := 0;
        end if;
        tauo := D / T;
        tauc := p1 * exp(p2 * tauo) + p3 * exp(p4 * tauo);
        kp := (tauo * (2 + tauo / 2) * (3 * tauc + tauo / 2) - 2 * tauc ^ 3 - 3 * tauc ^ 2 * tauo) / 2 / (tauc + tauo / 2) ^ 3;
        taui := (tauo * (2 + tauo / 2) * (3 * tauc + tauo / 2) - 2 * tauc ^ 3 - 3 * tauc ^ 2 * tauo) / 2 / tauo / (1 + tauo);
        taud := (3 * tauc ^ 2 * tauo + tauo ^ 2 / 2 * (3 * tauc + tauo / 2) - 2 * (1 + tauo) * tauc ^ 3) / (tauo * (2 + tauo / 2) * (3 * tauc + tauo / 2) - 2 * tauc ^ 3 - 3 * tauc ^ 2 * tauo);
        K := kp / mu;
        Ti := taui * T;
        Td := taud * T;
        N := kk / K * (1 + nn) - 1;
      end PID_VAP18NpropCinf;

      function PID_prop "prop rule for min pmarg=35 (default) or 25,30,40,45, ISA PID"
        input Real mu, T, D;
        input Integer pmargmin = Constants.pmminRp;
        output Real K, Ti, Td, N;
      protected
        Real theta, xqi, xpm, x, kappa, tau, bc0, bc1, bc2, ac1, ac2;
      algorithm
        theta := D / (D + T);
        xqi := 0.5515032 * theta - 1.0133502 * theta ^ 2 + 3.2132622 * theta ^ 3 - 1.8428728 * theta ^ 4;
        if pmargmin == 25 then
          xpm := 0.517 * theta + 0.0074 - 0.005 * exp(-theta / 0.065);
        elseif pmargmin == 30 then
          xpm := 0.542 * theta + 0.009 - 0.006 * exp(-theta / 0.062);
        elseif pmargmin == 35 then
          xpm := 0.531 * theta + 0.0247 - 0.0245 * exp(-theta / 0.06);
        elseif pmargmin == 40 then
          xpm := 0.523 * theta + 0.0409 - 0.04 * exp(-theta / 0.055);
        elseif pmargmin == 45 then
          xpm := 0.499 * theta + 0.0652 - 0.06 * exp(-theta / 0.05);
        else
          xpm := 0;
        end if;
        
        if Constants.enforce_pmminRp then
           x := max(xqi, xpm);
        else
           x := xqi;
        end if;
        
        kappa := x / (1 - x);
        tau := kappa * T;
        bc0 := 1 / mu;
        bc1 := (4 * D ^ 3 * T ^ 2 + 32 * tau * D ^ 2 * T ^ 2 + D ^ 4 * T + 8 * tau * D ^ 3 * T - 24 * tau ^ 2 * D ^ 2 * T - 32 * tau ^ 3 * D * T + 16 * tau ^ 4 * T + 16 * tau ^ 4 * D) / (4 * mu * D ^ 2 * T * (2 * T + D));
        bc2 := (D ^ 4 * T ^ 2 + 8 * tau * D ^ 3 * T ^ 2 + 24 * tau ^ 2 * D ^ 2 * T ^ 2 - 32 * tau ^ 3 * D * T ^ 2 + 16 * tau ^ 4 * T ^ 2 - 32 * tau ^ 3 * D ^ 2 * T + 16 * tau ^ 4 * D * T + 8 * tau ^ 4 * D ^ 2) / (4 * mu * D ^ 2 * T * (2 * T + D));
        ac1 := (D ^ 4 * T + 8 * tau * D ^ 3 * T + 24 * tau ^ 2 * D ^ 2 * T + 32 * tau ^ 3 * D * T - 16 * tau ^ 4 * T - 16 * tau ^ 4 * D) / (4 * D ^ 2 * T * (2 * T + D));
        ac2 := 2 * tau ^ 4 / (D * T);
        K := (ac1 * bc1 - ac2 * bc0) / ac1 ^ 2;
        Ti := (ac1 * bc1 - ac2 * bc0) / (ac1 * bc0);
        Td := (ac2 ^ 2 * bc0 - ac1 * ac2 * bc1 + ac1 ^ 2 * bc2) / (ac1 ^ 2 * bc1 - ac1 * ac2 * bc0);
        N := (ac2 ^ 2 * bc0 - ac1 * ac2 * bc1 + ac1 ^ 2 * bc2) / (ac1 * ac2 * bc1 - ac2 ^ 2 * bc0);
      end PID_prop;

      function ISAPID_pclass_ppar_trule "tune ISA PID for process class c and parameter value p, with rule r"
        input Integer c(min = 1, max = Constants.nClasses) "process class (1-Constants.nClasses)";
        input Real p "process class parameter";
        input Integer r(min = 1, max = Constants.nRules) "tuning rule (1-Constants.nRules)";
        output Real K, Ti, Td, N, mu, T, D;
      algorithm
        (mu, T, D) := Functions.AnalyticalMoA.MoA_P_all(c, p);
        if r == 1 then
          (K, Ti, Td, N) := Functions.TuningRules.PID_KS88IAEr(mu, T, D);
        elseif r == 2 then
          (K, Ti, Td, N) := Functions.TuningRules.PID_TS95dir(mu, T, D);
        elseif r == 3 then
          (K, Ti, Td, N) := Functions.TuningRules.PID_LC04std(mu, T, D);
        elseif r == 4 then
          (K, Ti, Td, N) := Functions.TuningRules.PID_LC04small(mu, T, D);
        elseif r == 5 then
          (K, Ti, Td, N) := Functions.TuningRules.PID_AR10r(mu, T, D);
        elseif r == 6 then
          (K, Ti, Td, N) := Functions.TuningRules.PID_AMIGO(mu, T, D);
        elseif r == 7 then
          (K, Ti, Td, N) := Functions.TuningRules.PI_SIMC(mu, T, D);
        elseif r == 8 then
          (K, Ti, Td, N) := Functions.TuningRules.PID_VAP18NVAP(mu, T, D);
        elseif r == 9 then
          (K, Ti, Td, N) := Functions.TuningRules.PID_VAP18NpropCinf(mu, T, D);
        else
          (K, Ti, Td, N) := Functions.TuningRules.PID_prop(mu, T, D);
        end if;
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
      end ISAPID_pclass_ppar_trule;
    end TuningRules;
  end Functions;

  package BenchmarkProcesses "processes for benchamrk tests"
    package BaseClasses
      model Bench_P
      equation
  
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
      end Bench_P;
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
    end BaseClasses;
  
    model Bench_Pnom "(1,1) Pade approx of FOPDT with normalised delay p"
    //A:matrix([0,1],[(2*p-2)/p,(p-2)/p]);
    //b:matrix([0],[1]);
    //c:matrix([(2-2*p)/p,-1]);
      extends BaseClasses.Bench_P;
      parameter Real p = 0.1;
      input Real u;
      output Real y;
    protected
      final parameter Real D = p/(1-p);
      Real[2] x(each start = 0, each fixed = true);
    equation
      der(x[1]) = x[2];
      der(x[2]) = (x[1]*(2*p-2))/p+(x[2]*(p-2))/p+u;
      y = (x[1]*(2-2*p))/p-x[2];
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end Bench_Pnom;
  
    model Bench_P1 "1/(1+s)^n, process class 1 in Astrom & Hagglund, 2000"
      extends BaseClasses.Bench_P;
      parameter Real p = 1;
      input Real u;
      output Real y;
    protected
      final parameter Integer n = integer(p);
      Real[n] x(each start = 0, each fixed = true);
    equation
      der(x[1]) + x[1] = u;
      for i in 2:n loop
        der(x[i]) + x[i] = x[i - 1];
      end for;
      y = x[n];
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end Bench_P1;

    model Bench_P2 "1/((1+s)(1+ps)(1+p^2s)(1+p^3s), process class 2 in Astrom & Hagglund, 2000"
      extends BaseClasses.Bench_P;
      parameter Real p = 1;
      input Real u;
      output Real y;
    protected
      Real[4] x(each start = 0, each fixed = true);
    equation
      der(x[1]) + x[1] = u;
      for i in 2:4 loop
        p ^ (i - 1) * der(x[i]) + x[i] = x[i - 1];
      end for;
      y = x[4];
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
    end Bench_P2;

    model Bench_P3 "(1-p*s)/(1+s)^3, process class 3 in Astrom & Hagglund, 2000"
      extends BaseClasses.Bench_P;
      parameter Real p = 1;
      input Real u;
      output Real y;
    protected
      Real[3] x(each start = 0, each fixed = true);
      Real y1;
    equation
      der(x[1]) + x[1] = (1 + p) * u "y1/u = (1-ps)/(1+s) = (1+p)/(1+s)-p";
      y1 = x[1] - p * u;
      der(x[2]) + x[2] = y1;
      der(x[3]) + x[3] = x[2];
      y = x[3];
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.02));
    end Bench_P3;

    model Bench_P4 "exp(-s)/(1+ps), process class 4 in Astrom & Hagglund, 2000"
      extends BaseClasses.Bench_P;
      parameter Real p = 1;
      input Real u;
      output Real y;
    protected
      Real x(start = 0, fixed = true);
    equation
      p * der(x) + x = delay(u, 1);
      y = x;
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.02));
    end Bench_P4;

    model Bench_P5 "exp(-s)/(1+ps)^2, process class 5 in Astrom & Hagglund, 2000"
      extends BaseClasses.Bench_P;
      parameter Real p = 1;
      input Real u;
      output Real y;
    protected
      Real x[2](each start = 0, each fixed = true);
    equation
      p * der(x[1]) + x[1] = delay(u, 1);
      p * der(x[2]) + x[2] = x[1];
      y = x[2];
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.02));
    end Bench_P5;

  
  end BenchmarkProcesses;

  package ModelsAndControllers "models for tuning and controllers"
    model FOPDT "mu exp(-sD)/(1+sT)"
      parameter Real mu(fixed = false);
      parameter Real T(fixed = false);
      parameter Real D(fixed = false);
      input Real u;
      output Real y(start = 0, fixed = true);
    equation
      y + T * der(y) = mu * delay(u, D);
    end FOPDT;

    model ISAPID "K(1+1/sTi+sTd/(1+sTd/N)), linear, continuous-time"
      parameter Real K;
      parameter Real Ti;
      parameter Real Td;
      parameter Real N;
      input Real u;
      output Real y;
    protected
      Real xi(start = 0, fixed = true);
      Real xd(start = 0, fixed = true);
      Real yd;
    equation
      der(xi) = K / Ti * u;
/* sKTd/(1+sTd/N) = KN-KN/(1+sTd/N) */
      if abs(Td) > 1e-12 then
        xd + Td / N * der(xd) = K * N * u;
        yd = K * N * u - xd;
      else
        xd + der(xd) = 0;
        yd = 0;
      end if;
      y = xi + yd + K * u;
    end ISAPID;
  end ModelsAndControllers;

  package Loops "loops composed of controler and benchmark process"
    package BaseClasses
      partial model Loop_ISAPID_Bench_P
        input Real r "set point";
        input Real d "load disturbance";
        output Real y "controlled variable";
        output Real u "control signal";
        //protected
        replaceable model Bench_P = BenchmarkProcesses.Bench_P1 constrainedby BenchmarkProcesses.BaseClasses.Bench_P;
        parameter Real p;
        parameter Real K(fixed = false);
        parameter Real Ti(fixed = false);
        parameter Real Td(fixed = false);
        parameter Real N(fixed = false);
        /* model parameters used only for computing boundary values, not for simuation */
        parameter Real mu(fixed = false);
        parameter Real T(fixed = false);
        parameter Real D(fixed = false);
        Bench_P P(p = p);
        ModelsAndControllers.ISAPID C(K = K, Ti = Ti, Td = Td, N = N);
      equation
        y = P.y;
        u = C.y;
        C.u = r - y;
        P.u = u + d;
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
      end Loop_ISAPID_Bench_P;
    end BaseClasses;
    
      model Loop_ISAPID_Bench_Pnom
      extends BaseClasses.Loop_ISAPID_Bench_P(redeclare model Bench_P = BenchmarkProcesses.Bench_Pnom);
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end Loop_ISAPID_Bench_Pnom;

    model Loop_ISAPID_Bench_P1
      extends BaseClasses.Loop_ISAPID_Bench_P(redeclare model Bench_P = BenchmarkProcesses.Bench_P1);
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end Loop_ISAPID_Bench_P1;

    model Loop_ISAPID_Bench_P2
      extends BaseClasses.Loop_ISAPID_Bench_P(redeclare model Bench_P = BenchmarkProcesses.Bench_P2);
      annotation(
        experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
    end Loop_ISAPID_Bench_P2;

    model Loop_ISAPID_Bench_P3
      extends BaseClasses.Loop_ISAPID_Bench_P(redeclare model Bench_P = BenchmarkProcesses.Bench_P3);
    end Loop_ISAPID_Bench_P3;

    model Loop_ISAPID_Bench_P4
      extends BaseClasses.Loop_ISAPID_Bench_P(redeclare model Bench_P = BenchmarkProcesses.Bench_P4);
    end Loop_ISAPID_Bench_P4;

    model Loop_ISAPID_Bench_P5
      extends BaseClasses.Loop_ISAPID_Bench_P(redeclare model Bench_P = BenchmarkProcesses.Bench_P5);
    end Loop_ISAPID_Bench_P5;
  end Loops;

  package Test_AnalyticalMoA "models to test the FOPDT parametrisation via MoA"
    package BaseClasses
      partial model Test_MoA_Bench_P
        parameter Real[:] p = {1};
        Real yp[n], ym[n];
        //protected
        final parameter Integer n = size(p, 1);
        replaceable model Process = BenchmarkProcesses.Bench_P1 constrainedby BenchmarkProcesses.BaseClasses.Bench_P;
        replaceable function MoA = Functions.AnalyticalMoA.MoA_P1;
        Process[n] P(p = p);
        ModelsAndControllers.FOPDT[n] M;
      initial equation
        for i in 1:n loop
          (M[i].mu, M[i].T, M[i].D) = MoA(p[i]);
        end for;
      equation
        for i in 1:n loop
          yp[i] = P[i].y;
          ym[i] = M[i].y;
          P[i].u = if time < 0.1 then 0 else 1;
          M[i].u = if time < 0.1 then 0 else 1;
        end for;
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
          experiment(StartTime = 0, StopTime = 30, Tolerance = 1e-06, Interval = 0.06));
      end Test_MoA_Bench_P;
    end BaseClasses;

    model Test_MoA_Bench_P1
      extends BaseClasses.Test_MoA_Bench_P(redeclare model Process = BenchmarkProcesses.Bench_P1, redeclare function MoA = Functions.AnalyticalMoA.MoA_P1, p = {1, 2, 4, 8, 10});
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        experiment(StartTime = 0, StopTime = 30, Tolerance = 1e-06, Interval = 0.06));
    end Test_MoA_Bench_P1;

    model Test_MoA_Bench_P2
      extends BaseClasses.Test_MoA_Bench_P(redeclare model Process = BenchmarkProcesses.Bench_P2, redeclare function MoA = Functions.AnalyticalMoA.MoA_P2, p = {0.1, 0.2, 0.4, 0.8, 1.0});
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        experiment(StartTime = 0, StopTime = 30, Tolerance = 1e-06, Interval = 0.06));
    end Test_MoA_Bench_P2;

    model Test_MoA_Bench_P3
      extends BaseClasses.Test_MoA_Bench_P(redeclare model Process = BenchmarkProcesses.Bench_P3, redeclare function MoA = Functions.AnalyticalMoA.MoA_P3, p = {0.2, 0.5, 1, 2, 3});
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        experiment(StartTime = 0, StopTime = 30, Tolerance = 1e-06, Interval = 0.06));
    end Test_MoA_Bench_P3;

    model Test_MoA_Bench_P4
      extends BaseClasses.Test_MoA_Bench_P(redeclare model Process = BenchmarkProcesses.Bench_P4, redeclare function MoA = Functions.AnalyticalMoA.MoA_P4, p = {0.1, 0.5, 1, 5, 10});
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        experiment(StartTime = 0, StopTime = 30, Tolerance = 1e-06, Interval = 0.06));
    end Test_MoA_Bench_P4;

    model Test_MoA_Bench_P5
      extends BaseClasses.Test_MoA_Bench_P(redeclare model Process = BenchmarkProcesses.Bench_P5, redeclare function MoA = Functions.AnalyticalMoA.MoA_P5, p = {0.1, 0.5, 1, 5, 10});
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})),
        experiment(StartTime = 0, StopTime = 80, Tolerance = 1e-06, Interval = 0.060015));
    end Test_MoA_Bench_P5;
  end Test_AnalyticalMoA;

  package BenchmarkTests "benchmark tests: closed-loop responses to unit load disturbance step at t=0"
    package BaseClasses
      partial model LoadDistRes_Loop_P
      protected
        replaceable model Bench_Loop = Loops.Loop_ISAPID_Bench_P1 constrainedby Loops.BaseClasses.Loop_ISAPID_Bench_P;
        Bench_Loop[np, nr] L(p = Functions.Internal.vec2mat_byrow(ppars, nr));
        parameter Integer pclass = 1;
        parameter Integer[:] rules = 1:Constants.nRules;
        parameter Real[:] ppars = {2, 4, 8, 10};
        parameter Real csvSampleTime = 0.1;
        final parameter Integer np = size(ppars, 1);
        final parameter Integer nr = size(rules, 1);
        Real[np, nr] ISE(each start = 0, each fixed = true);
        Real[np, nr] IAE(each start = 0, each fixed = true);
        Real[np, nr] ITAE(each start = 0, each fixed = true);
        Real[np, nr] LRPI(each start = 0, each fixed = true);
        /* Visioli 2012, Ad=1 */
        String[np] csvName_resp;
        String csvName_TFparams;
        String csvLine;
      public
        Real y[np, nr];
      initial equation
        for i in 1:np loop
          for j in 1:nr loop
            (L[i, j].K, L[i, j].Ti, L[i, j].Td, L[i, j].N, L[i, j].mu, L[i, j].T, L[i, j].D) = Functions.TuningRules.ISAPID_pclass_ppar_trule(pclass, ppars[i], rules[j]);
          end for;
        end for;
      equation
        for i in 1:np loop
          for j in 1:nr loop
            L[i, j].r = 0;
            L[i, j].d = if time < Modelica.Constants.small then 0 else 1;
            y[i, j] = L[i, j].y;
            der(ISE[i, j]) = L[i, j].y ^ 2;
            der(IAE[i, j]) = abs(L[i, j].y);
            der(ITAE[i, j]) = time * abs(L[i, j].y);
          end for;
        end for;
      algorithm
        when initial() then
          csvLine := "time";
          for i in 1:nr loop
            csvLine := csvLine + ",y" + String(i, "02d");
          end for;
          for i in 1:np loop
            csvName_resp[i] := "class_" + String(pclass, "02d") + "_param_" + String(ppars[i]) + ".csv";
            Modelica.Utilities.Files.remove(csvName_resp[i]);
            Modelica.Utilities.Streams.print(csvLine, csvName_resp[i]);
          end for;
        end when;
        when sample(0, csvSampleTime) then
          for i in 1:np loop
            csvLine := String(time);
            for j in 1:nr loop
              csvLine := csvLine + "," + String(y[i, j]);
            end for;
            Modelica.Utilities.Streams.print(csvLine, csvName_resp[i]);
          end for;
        end when;
        when terminal() then
          for i in 1:np loop
            Modelica.Utilities.Streams.close(csvName_resp[i]);
          end for;
          Functions.Internal.write_matrix_to_csv(ISE, "ISE_pclass_" + String(pclass, "02d") + ".csv");
          Functions.Internal.write_matrix_to_csv(IAE, "IAE_pclass_" + String(pclass, "02d") + ".csv");
          Functions.Internal.write_matrix_to_csv(ITAE, "ITAE_pclass_" + String(pclass, "02d") + ".csv");
          for i in 1:np loop
            for j in 1:nr loop
              LRPI[i, j] := 27 * L[i, j].mu * L[i, j].D ^ 2 / 4 / (2 * L[i, j].T + L[i, j].D) / IAE[i, j];
            end for;
          end for;
          Functions.Internal.write_matrix_to_csv(LRPI, "LRPI_pclass_" + String(pclass, "02d") + ".csv");
          csvName_TFparams := "class_" + String(pclass, "02d") + "_TFparams.csv";
          Modelica.Utilities.Files.remove(csvName_TFparams);
          csvLine := "param,rule,K,Ti,Td,N,mu,T,D";
          Modelica.Utilities.Streams.print(csvLine, csvName_TFparams);
          for i in 1:np loop
            for j in 1:nr loop
              csvLine := String(ppars[i]) + "," + String(rules[j]) + "," + String(L[i, j].K) + "," + String(L[i, j].Ti) + "," + String(L[i, j].Td) + "," + String(L[i, j].N) + "," + String(L[i, j].mu) + "," + String(L[i, j].T) + "," + String(L[i, j].D);
              Modelica.Utilities.Streams.print(csvLine, csvName_TFparams);
            end for;
          end for;
          Modelica.Utilities.Streams.close(csvName_TFparams);
        end when;
        annotation(
          Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
      end LoadDistRes_Loop_P;
      annotation(
        Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
    end BaseClasses;
    
    model LoadDistRes_Loop_Pnom
      extends BaseClasses.LoadDistRes_Loop_P(redeclare model Bench_Loop = Loops.Loop_ISAPID_Bench_Pnom, pclass = 0, ppars = {0.05,0.1,0.2,0.4,0.6,0.9}, csvSampleTime = 0.05, rules = 1:Constants.nRules);
      annotation(
        experiment(StartTime = 0, StopTime = 30, Tolerance = 1e-06, Interval = 0.02),
        __OpenModelica_simulationFlags(s = "dassl"));
    end LoadDistRes_Loop_Pnom;

    model LoadDistRes_Loop_P1
      extends BaseClasses.LoadDistRes_Loop_P(redeclare model Bench_Loop = Loops.Loop_ISAPID_Bench_P1, pclass = 1, ppars = {2, 3, 4, 6, 8, 10}, csvSampleTime = 0.05, rules = 1:Constants.nRules);
      annotation(
        experiment(StartTime = 0, StopTime = 60, Tolerance = 1e-06, Interval = 0.02),
        __OpenModelica_simulationFlags(s = "dassl"));
    end LoadDistRes_Loop_P1;

    model LoadDistRes_Loop_P2 "Use preferably the rungeKuttaSsc solver"
      extends BaseClasses.LoadDistRes_Loop_P(redeclare model Bench_Loop = Loops.Loop_ISAPID_Bench_P2, pclass = 2, ppars = {0.1, 0.2, 0.4, 0.6, 0.8, 0.9}, csvSampleTime = 0.05, rules = 1:Constants.nRules);
      annotation(
        experiment(StartTime = 0, StopTime = 20, Tolerance = 1e-06, Interval = 0.01),
        __OpenModelica_simulationFlags(s = "dassl"));
    end LoadDistRes_Loop_P2;

    model LoadDistRes_Loop_P3 "Use preferably the rungeKuttaSsc solver"
      extends BaseClasses.LoadDistRes_Loop_P(redeclare model Bench_Loop = Loops.Loop_ISAPID_Bench_P3, pclass = 3, ppars = {0.1, 0.2, 0.5, 1, 2, 3}, csvSampleTime = 0.05, rules = 1:Constants.nRules);
      annotation(
        experiment(StartTime = 0, StopTime = 30, Tolerance = 1e-06, Interval = 0.02),
        __OpenModelica_simulationFlags(s = "dassl"));
    end LoadDistRes_Loop_P3;

    model LoadDistRes_Loop_P4
      extends BaseClasses.LoadDistRes_Loop_P(redeclare model Bench_Loop = Loops.Loop_ISAPID_Bench_P4, pclass = 4, ppars = {0.1, 0.2, 0.5, 2, 5, 10}, csvSampleTime = 0.05, rules = 1:Constants.nRules);
      annotation(
        experiment(StartTime = 0, StopTime = 30, Tolerance = 0.0001, Interval = 0.12),
        __OpenModelica_simulationFlags(s = "dassl"));
    end LoadDistRes_Loop_P4;

    model LoadDistRes_Loop_P5
      extends BaseClasses.LoadDistRes_Loop_P(redeclare model Bench_Loop = Loops.Loop_ISAPID_Bench_P5, pclass = 5, ppars = {0.1, 0.2, 0.5, 2, 5, 10}, csvSampleTime = 0.05, rules = 1:Constants.nRules);
      annotation(
        experiment(StartTime = 0, StopTime = 79.8, Tolerance = 0.0001, Interval = 0.149346),
        __OpenModelica_simulationFlags(s = "dassl"));
    end LoadDistRes_Loop_P5;
  end BenchmarkTests;
  annotation(
    Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})));
end TuneForLD_Qshape;