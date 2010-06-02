use v6;
use Test;

plan *;

BEGIN { @*INC.push: 'lib' };
use Math::RungeKutta;

# simple differential equation:
# parameter: x
# value:     y
# y = x**2
# dy/dy = 2x
# y0 = 0

my $last = 0;
sub record-last($t, @y) { $last = @y[0] };

lives_ok { rk-integrate(
    :from(0),
    :to(3),
    :initial[0],
    :derivative(-> $x, @y { 2 * $x }),
    :do(&record-last),
    :order(4)
) }, 'lives through a simple RK4 integration';

is_approx($last, 3**2, 'and produced a good approximation of x**2 with x = 3');

done_testing;
