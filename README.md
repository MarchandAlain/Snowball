# Snowball
I present here a novel snowball algorithm to compute the minimal enclosing disk for a set of 2-D points. At each step, the set supporting a minimal disk for a sample of four points is determined, and the first one or two points not included in the disk are substituted to the unnecessary points in the sample, always searching forward (cyclically) in the set. I demonstrate that the algorithm always yields the minimal enclosing disks for the set. It runs in linear time and requires little more than two inclusion tests (distance calculations) per point. Extensions to higher dimensions and ellipsoids may be considered.