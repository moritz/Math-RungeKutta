module Math::RungeKutta;

# the actual integation loop
sub rk-integrate (
    :$from,
    :$to,
    :@initial,
    :$step = ($from - $to) / 100,
    :&derivative,
    :&do = sub ($parameter, @values) { say "$parameter: @values.join(', ')" },
    :$order = 4,
) is export {
    my $parameter = $from;
    my @values = @initial;
    while $parameter < $to {
        @values = runge-kutta(
                :$order,
                :@values,
                :$parameter,
                :&derivative,
                :$step,
        );
        $parameter = @values.shift;
        do($parameter, @values);
    }
}

sub runge-kutta(:@values, :&derivative,
                :$parameter, :$step = 0.01, :$order = 4) is export {
    rk(:values(@values), :&derivative, :$parameter, :$step, :$order);
}

# Euler's algorithm
multi sub rk(:@values, :&derivative, :$parameter, :$step = 0.1, :$order where 1) {
    my @new_y = @values >>+<< ($step <<*<< derivative($parameter, @values));
    return ($parameter + $step, @new_y);
}

# Heun's method
multi sub rk(:@values, :&derivative, :$parameter, :$step = 0.1, :$order where 2) {
    my @k1 = $step <<*<< derivative($parameter, @values);
    my @k2 = $step <<*<< derivative($parameter+$step, @values >>+<< @k1);
    my @result = @values >>+<< (1/2) <<*<< (@k1 >>+<< @k2);
    return ($parameter + $step, @result);
}

# classical RK4
multi sub rk(:@values, :&derivative, :$parameter, :$step = 0.1, :$order where 4) {
    my @k1 = $step <<*<< derivative($parameter, @values);
    my @k2 = $step <<*<< derivative($parameter + $step/2, [@values >>+<< @k1 >>/>> 2]);
    my @k3 = $step <<*<< derivative($parameter + $step/2, [@values >>+<< @k2 >>/>>2]);
    my @k4 = $step <<*<< derivative($parameter + $step,   [@values >>+<< @k3]);
    my @result = @values >>+<< ((1/6) <<*<< (@k1 >>+<< @k4))
                         >>+<< ((1/3) <<*<< (@k2 >>+<< @k3));
    return ($parameter + $step, @result);
}

# vim: ft=perl6
