# Introduction

Runge-Kutta integration for systems of ordinary, linear differential
equations.

WARNING: this module is not yet thoroughly tested. Use it at your own risk.
Bug reports, tests and patches are welcome!

Let's say you have a differential equation for the function f(t),
with the equation

    df/dt = f(t)^2 + t

and the initial value f(t=0) = 1;

In the context of this module, we call df/dt the "derivative",
t the "parameter"

You'd solve that numerically with this Perl 6 code:

    # synopsis.pl
    use Math::RungeKutta;
    
    # function that calculates the derivative that
    # Math::RungeKutta will integrate
    sub d($t, @values) { @values[0]**2 + $t}
    
    # that's a function that gets called with the
    # current values after each integration step
    sub callback($t, @values) { say "$t @values[0]" };
    
    my @initial = 1;
    
    adaptive-rk-integrate(
        :from(0),
        :to(0.5),
        :@initial,
        :derivative(&d),
        :do(&callback)
    );

And then look at the result:

    $ PERL6LIB=lib perl6 synopsis.pl | xmgrace -nxy -

The interfaces is inspired by Perl 5 module Math::RungeKutta, to be found at
<http://search.cpan.org/perldoc?Math::RungeKutta>.

# Build Status

[![Build Status](https://travis-ci.org/moritz/Math-RungeKutta.svg)](https://travis-ci.org/moritz/Math-RungeKutta)

# License

This module is licensed under the [Artistic License version 2.0](https://opensource.org/licenses/Artistic-2.0).


Its accompanying tests and examples are public domain, as defined by the [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/).
