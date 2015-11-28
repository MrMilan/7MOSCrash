package MOS7
  model LC
    Modelica.Electrical.Analog.Basic.Inductor inductor1(L = 0.1, i(fixed = true)) annotation(Placement(visible = true, transformation(origin = {40, 2}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Resistor resistor1(R = 1) annotation(Placement(visible = true, transformation(origin = {6, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Capacitor capacitor1(C = 0.001, v(start = 1, fixed = true)) annotation(Placement(visible = true, transformation(origin = {-14, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Electrical.Analog.Basic.Ground ground1 annotation(Placement(visible = true, transformation(origin = {30, -36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(resistor1.n, ground1.p) annotation(Line(points = {{16, -26}, {30, -26}}, color = {0, 0, 255}));
    connect(ground1.p, inductor1.n) annotation(Line(points = {{30, -26}, {40, -26}, {40, -8}}, color = {0, 0, 255}));
    connect(resistor1.p, capacitor1.p) annotation(Line(points = {{-4, -26}, {-14, -26}, {-14, -8}, {-14, -8}}, color = {0, 0, 255}));
    connect(capacitor1.n, inductor1.p) annotation(Line(points = {{-14, 12}, {-14, 22}, {40, 22}, {40, 12}}, color = {0, 0, 255}));
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end LC;

  model advectionDiscretized
    // u_t + u_x = 0
    parameter Real L = 1;
    constant Integer N = 1000;
    parameter Real dx = L / (N - 1);
    parameter Real[N] x = array(i * dx for i in 0:N - 1);
    Real[N] u(each start = 0);
    parameter Real c = 1;
  initial equation

  equation
    for i in 2:N - 1 loop
      // the equation:
      der(u[i]) + c * (u[i + 1] - u[i - 1]) / (2 * dx) = 0;
    end for;
    u[1] = sin(2 * 3.14159 * time);
    //left BC
    u[N] = 2 * u[N - 1] - u[N - 2];
    //extrapolation in the last node
    annotation(experiment(Interval = 0.002));
  end advectionDiscretized;

  model ode
    Real x(start = 1, fixed = true);
    parameter Real p = 3;
  equation
    der(x) = -p * x;
  end ode;

  model koreny
    Real x, y, z;
    parameter Real ep = 0.005;
    Boolean b;
  equation
    x = time;
    y = 68 - 40 * x - x ^ 2 + 2 * x ^ 3;
    if x > (-ep) then
      z = sin(2 * ep) + log(ep) + ep + 8;
      b = true;
    else
      z = sin(2 * x) + log(-x) + x + 8;
      b = false;
    end if;
    annotation(experiment(StartTime = -10, StopTime = 10, Tolerance = 1e-06, Interval = 0.04));
  end koreny;

  model inicializace1
    Real x(start = 0);
  initial equation
    0 = 68 - 40 * x - x ^ 2 + 2 * x ^ 3;
  equation
    der(x) = -x;
  end inicializace1;

  model Rc
    Modelica.Electrical.Analog.Basic.Ground ground1 annotation(Placement(visible = true, transformation(origin = {-16, -38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor resistor1(R = 1) annotation(Placement(visible = true, transformation(origin = {-56, 2}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Capacitor capacitor1(C = 1, v(start = 10, fixed = true)) annotation(Placement(visible = true, transformation(origin = {34, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  equation
    connect(capacitor1.p, ground1.p) annotation(Line(points = {{34, -12}, {34, -12}, {34, -28}, {-16, -28}, {-16, -28}}, color = {0, 0, 255}));
    connect(resistor1.n, ground1.p) annotation(Line(points = {{-56, -8}, {-54, -8}, {-54, -28}, {-16, -28}, {-16, -28}}, color = {0, 0, 255}));
    connect(resistor1.p, capacitor1.n) annotation(Line(points = {{-56, 12}, {-56, 12}, {-56, 40}, {34, 40}, {34, 8}, {34, 8}}, color = {0, 0, 255}));
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  end Rc;

  model crash
    parameter Real m(start = 76, fixed = true);
    parameter Real v(start = 65 / 3.6, fixed = false);
    Real a(start = -30, fixed = false);
    Real E, s;
  initial equation
    der(v) = 0;
  equation
    E = m * v ^ 2 / 2;
    v = s / time;
    der(v) = a;
  end crash;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), uses(Modelica(version = "2.2.2")));
end MOS7;