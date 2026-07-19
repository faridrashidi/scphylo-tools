from types import SimpleNamespace


class _FakeGurobiExpression:
    """Minimal symbolic expression used by optional-solver adapter tests."""

    def __init__(self, value=0):
        self.X = value

    def _expression(self, other=None):
        return _FakeGurobiExpression()

    __add__ = _expression
    __eq__ = _expression
    __ge__ = _expression
    __le__ = _expression
    __mul__ = _expression
    __radd__ = _expression
    __rsub__ = _expression
    __rmul__ = _expression
    __sub__ = _expression

    def __neg__(self):
        return _FakeGurobiExpression()


def fake_gurobi(value_for_variable=lambda model, index, name: 0):
    """Return a tiny Gurobi-compatible API for wrapper tests.

    The fake deliberately does not optimize anything. It records the model calls and
    gives each variable a caller-selected solution value so tests can exercise the
    Python adapter deterministically without Gurobi or a solver license.
    """
    models = []

    class Model:
        def __init__(self, name):
            self.name = name
            self.Params = SimpleNamespace()
            self.variables = []
            self.constraints = []
            self.objective = None
            models.append(self)

        def addVar(self, **kwargs):
            index = len(self.variables)
            variable = _FakeGurobiExpression(
                value_for_variable(self.name, index, kwargs.get("name"))
            )
            self.variables.append(variable)
            return variable

        def addConstr(self, constraint):
            self.constraints.append(constraint)

        def update(self):
            return None

        def setObjective(self, objective, sense=None):
            self.objective = (objective, sense)

        def optimize(self):
            return None

    api = SimpleNamespace(
        GRB=SimpleNamespace(BINARY="binary", MAXIMIZE="maximize", MINIMIZE="minimize"),
        Model=Model,
        quicksum=lambda values: sum(values, _FakeGurobiExpression()),
    )
    api.models = models
    return api
