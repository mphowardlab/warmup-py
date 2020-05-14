r"""Concepts of molecular simulation in Python.

Two interfaces for classes are outlined but are incomplete:

.. autosummary::
    :nosignatures:

    warmup.LennardJones
    warmup.RadialDistribution

Unit tests are provided for both classes. You can run them from the project
directory using::

    python3 -m unittest

You can also run specific tests one at a time::

    python3 -m unittest test.test_LennardJones

Initially, all tests will fail! You will need to implement the details of these
classes in order to get the tests to pass.

Start with the :py:class:`LennardJones` pair potential and implement the
:py:meth:`~LennardJones.energy` and :py:meth:`~LennardJones.force` methods.
Make sure you pay attention to what types of inputs are allowed or returned,
and special cases of return values. Then, you can move onto
:py:class:`RadialDistribution`, which is a more involved class used for
analysis.

"""

class LennardJones:
    r"""Lennard-Jones pair potential.

    The prototypical pairwise interaction consisting of a steep repulsive core
    and an attractive tail describing dispersion forces. The functional form
    of the potential is:

    .. math::

        u(r) = \begin{cases}
               4 \varepsilon\left[\left(\dfrac{\sigma}{r}\right)^{12}
               - \left(\dfrac{\sigma}{r}\right)^6 \right], & r \le r_{\rm cut} \\
               0, & r > r_{\rm cut}
               \end{cases}

    where :math:`r` is the distance between the centers of two particles,
    :math:`\varepsilon` sets the strength of the attraction, and
    :math:`\sigma` sets the length scale of the interaction. (Typically,
    :math:`\sigma` can be regarded as a particle diameter.) The potential
    is truncated to zero at :math:`r_{\rm cut}` for computational efficiency,
    and a good accuracy for thermodynamic properties is usually achieved when
    :math:`r_{\rm cut} \ge 3\sigma`.

    In molecular dynamics (MD) simulations, the forces on the particles are what
    is actually required. Forces are computed from the gradient of :math:`u(r)`:

    .. math::

        \mathbf{F}(r) = \begin{cases}
                        -\nabla u(r), & r \le r_{\rm cut} \\
                        0, & r > r_{\rm cut}
                        \end{cases}

    A similar truncation scheme is applied to :math:`\mathbf{F}` as for
    :math:`u`. This implies that the energy should be shifted to zero at
    :math:`r_{\rm cut}` by subtracting :math:`u(r_{\rm cut})`. However, this
    distinction is often not made in MD simulations unless thermodynamic
    properties based on the energy are being computed. Caution must be taken
    if MD results with this scheme are compared to Monte Carlo (MC) results,
    which are sensitive to whether :math:`u` is shifted or not.

    If the Lennard-Jones potential is truncated and shifted at its minimum
    :math:`r_{\rm cut} = 2^{1/6}\sigma`, the interactions are purely repulsive.
    (The forces are always positive, :math:`|\mathbf{F}| \ge 0`.)
    This special case is often used as an approximation of nearly hard spheres,
    where it is referred to as the Weeks--Chandler--Andersen potential based
    on its role in their perturbation theory of the liquid state.

    Args:
        epsilon (float): Interaction energy.
        sigma (float): Interaction length.
        rcut (float): Truncation distance.

    Attributes:
        epsilon (float): Interaction energy.
        sigma (float): Interaction length.
        rcut (float): Truncation distance.

    """
    def __init__(self, epsilon, sigma, rcut):
        pass

    def energy(self, r):
        r"""Evaluate potential energy.

        Args:
            r (float or array_like): Pair distances.

        Returns:
            float or array_like: Energy at the pair distances.

            If ``r`` is 0, the energy is :py:obj:`numpy.inf`.

        """
        pass

    def force(self, r):
        r"""Evaluate force.

        Args:
            r (float or array_like): Pair distances

        Returns:
            float or array_like: Magnitude of force at the pair distances.

            If ``r`` is 0, the force is :py:obj:`numpy.inf`.

        """
        pass

class RadialDistribution:
    r"""Radial distribution function calculator.

    The radial distribution function measures the pairwise distance correlations
    between particles. (For this reason, it is also called a pair correlation
    function.) In an ideal gas, particles are noninteracting and randomly
    distributed. However, in real fluids or solids, there tends to be structural
    correlations between particles. For example, two hard particles cannot
    overlap each other, so they must always be separated by at least the sum
    of their radii.

    Suppose that we have :math:`N` indistinguishable particles in the system.
    These particles are interacting by a total potential energy function
    :math:`U` that depends on all of their coordinates :math:`\mathbf{r}^N`.
    The probability of observing two particles at certain positions
    :math:`\mathbf{r}_1` and :math:`\mathbf{r}_2` (with permutation of the
    tags), which we call the two-particle correlation function, is:

    .. math::

        \rho^{(2)}(\mathbf{r}_1,\mathbf{r}_2)
            = \frac{N(N-1)}{Z} \int {\rm d}\mathbf{r}^{N-2} e^{-\beta U(\mathbf{r}^N)}

    where :math:`Z` is a partition function. We will assume that the system is
    *homogeneous* so that the correlations depend on distances between particles
    rather than absolute coordinates, :math:`\mathbf{r} = \mathbf{r}_2 -
    \mathbf{r}_1`. We can define the pair correlation function as an ensemble
    average:

    .. math::

        g(\mathbf{r}) = \rho^{-1} \left\langle
            \sum_{i=2}^N \delta(\mathbf{r}-\mathbf{r}_i)\right\rangle

    where :math:`\rho = N/V` is the *average* density of :math:`N` particles in
    volume :math:`V`. The sum of delta functions effectively counts how many
    particles are at a certain separation from particle 1.

    .. note::
        In an ideal gas, this average can be taken exactly because :math:`U=0`
        and :math:`g_{\rm ig}(\mathbf{r}) = 1 - 1/N`. Another consequence of
        the definition of :math:`g(\mathbf{r})` is that:

        .. math::

            \int {\rm d}\mathbf{r} \rho g(\mathbf{r}) = N-1

    We will further assume that the correlations are isotropic so that they
    only depend on :math:`r = |\mathbf{r}|`. We can then compute the radial
    distribution function :math:`g(r)` by counting the number of particles
    :math:`\Delta N(r)` within a shell of volume :math:`\Delta V`:

    .. math::

        g(r) = \frac{\langle \Delta N(r) \rangle}{\rho \Delta V(r)}

    The bins are usually taken to be spherical shells of radial width
    :math:`\Delta r`. To improve the averaging of :math:`\Delta N`, multiple
    configurations are usually accumulated into an average, and each particle
    within the configuration is used as an origin for the calculation. However,
    only unique pairs of particles are counted because, e.g., the distance from
    particle 2 to 1 does not contribute new information if particle 1 to 2 has
    already been counted.

    When computing the distances between two particles, most molecular
    simulations employ *periodic boundary conditions* (as in the PAC-MAN arcade
    game!). In this convention, particles that exit one edge of the simulation
    volume reenter from the opposite side. We will focus on the one-dimensional
    problem, with the coordinate :math:`0 \le x < L`. In Python pseudocode,
    the particle position should be *wrapped* as::

        x -= L*np.floor(x/L)

    The periodic boundaries are equivalent to saying that the simulation volume
    is surrounded by copies of itself, and images of a particle move identically
    to the particle in the simulation volume. Since this system extends in
    principle to infinity, pairs of particles are then defined to interact only
    with their *nearest image*. In Python pseudocode, the vector between two
    particles should be wrapped as::

        dx -= L*np.round(dx/L)

    .. warning::

        Within this convention, it is not recommended to include pair distances
        :math:`r > L/2`. This calculator should raise an error if this would
        occur.

    Args:
        rmax (float): Maximum distance to compute distribution.
        dr (float): Bin width for shells.
        rmin (float, optional): Minimum distance to compute distribution.
            Defaults to 0.0.

    Attributes:
        rmax (float): Maximum distance to compute distribution.
        dr (float): Bin width for shells.
        rmin (float): Minimum distance to compute distribution.
        edges (:py:obj:`numpy.ndarray`): Edges of the bins for the distribution.
        centers (:py:obj:`numpy.ndarray`): Centers of the bins for the distribution.

    If ``dr`` does not evenly divide the range from ``rmin`` to ``rmax``,
    it will be rounded to the closest compatible value.

    """
    def __init__(self, rmax, dr, rmin=0.0):
        pass

    def accumulate(self, L, points):
        r"""Accumulate a new entry in the average.

        Args:
            L (float): Edge length of the cubic simulation volume.
            points (array_like): An ``Nx3`` array of points.

        Raises:
            ValueError: If ``rmax >= L/2``.

        The sample is accumulated into the radial distribution average.

        """
        pass

    def reset(self):
        r"""Reset the accumulator.

        The sample accumulator is reset to zeros.

        """
        pass

    @property
    def rdf(self):
        r""":py:obj:`numpy.ndarray`: Average radial distribution function."""
        pass
