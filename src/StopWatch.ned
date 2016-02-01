simple StopWatch
{
	parameters:
		double tti @unit(s);
		double epsilon @unit(s);
		double runtime @unit(s);
		double initOffset @unit(s);
        @display("p=16,16;i=device/clock");
}
