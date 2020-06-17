#define PY_SSIZE_T_CLEAN
#include <Python.h>

typedef struct {
  PyObject_HEAD
} MassActionObject;

/**
 * Choose function implemented as specified by the following StackOverflow post:
 * https://stackoverflow.com/questions/24294192/computing-the-binomial-coefficient-in-c
 */
static int choose(int n, int k)
{
  int ans = 1;
  k = k > n - k ? n - k : k;

  for (int j = 1; j <= k; ++j, --n) {
    if (n % j == 0) {
      ans *= n / j;
    } else {
      if (ans % j == 0) {
        ans = ans / j * n;
      } else {
        ans = (ans * n) / j;
      }
    }
  }

  return ans;
}

static PyObject *
MassAction_propensity(PyObject *self, PyObject *args)
{
  PyObject *state;
  if (!PyArg_ParseTuple(args, "O", &state)) {
    return NULL;
  }

  // Import 'stocal.structures.multiset'
  PyObject *structures = PyImport_ImportModule("stocal.structures");
  if (structures == NULL) {
    return NULL;
  }
  PyObject *multiset = PyObject_GetAttrString(structures, "multiset");

  // Ensure 'state' is a 'multiset'
  if (!PyObject_IsInstance(state, multiset)) {
    PyObject *fnargs = Py_BuildValue("(O)", state);
    PyObject *temp = PyObject_CallObject(multiset, fnargs);
    Py_DECREF(fnargs);

    if (temp == NULL) {
      return NULL;
    }

    Py_DECREF(state);
    state = temp;
  }

  // Check if state is empty
  if (PyDict_Size(state) == 0) {
    Py_DECREF(state);
    return PyFloat_FromDouble(0.0);
  }

  // Calculate propensity
  PyObject *constant = PyObject_GetAttrString(self, "constant");
  if (constant == NULL) {
    PyErr_SetString(PyExc_Exception, "constant not defined");
    Py_DECREF(state);
    return NULL;
  }

  // Define infinity
  PyObject *infstr = PyUnicode_FromString("inf");
  PyObject *inf = PyFloat_FromString(infstr);
  Py_DECREF(infstr);

  // Check if constant was infinity
  if (PyObject_RichCompareBool(constant, inf, Py_EQ)) {
    Py_DECREF(state);
    return inf;
  }

  PyObject *reactants = PyObject_GetAttrString(self, "reactants");
  if (reactants == NULL) {
    PyErr_SetString(PyExc_Exception, "reactants not defined");
    Py_DECREF(state);
    return NULL;
  }

  double a = PyFloat_AsDouble(constant);
  PyObject *s, *n, *stateval;
  Py_ssize_t pos = 0;

  while (PyDict_Next(reactants, &pos, &s, &n)) {
    long reactlng = PyLong_AsLong(n);
    stateval = PyDict_GetItem(state, s);

    long statelng = 0;
    if (stateval != NULL) {
      statelng = PyLong_AsLong(stateval);
    }

    /*
    if (stateval == NULL) {
      const char *sstr = PyUnicode_AsUTF8(s);
      PyErr_Format(PyExc_Exception, "key '%s' not found in state", sstr);
      return NULL;
    }
    */

    a *= choose(statelng, reactlng);
  }

  Py_DECREF(state);
  return PyFloat_FromDouble(a);
}

static PyMethodDef MassActionMethods[] = {
  {"propensity", (PyCFunction) MassAction_propensity, METH_VARARGS, "Reaction propensity for the given state."},
  {NULL} // Sentinel
};

static PyTypeObject MassActionType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "cstocal.MassAction",
  .tp_doc = "MassAction C implementation.",
  .tp_basicsize = sizeof(MassActionObject),
  .tp_itemsize = 0,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_new = PyType_GenericNew,
  .tp_methods = MassActionMethods
};

static PyModuleDef cstocalmodule = {
  PyModuleDef_HEAD_INIT,
  .m_name = "cstocal",
  .m_doc = "Stocal extension module.",
  .m_size = -1
};

PyMODINIT_FUNC
PyInit_cstocal(void)
{
  if (PyType_Ready(&MassActionType) < 0) {
    return NULL;
  }

  PyObject *mod = PyModule_Create(&cstocalmodule);
  if (mod == NULL) {
    return NULL;
  }

  Py_INCREF(&MassActionType);
  if (PyModule_AddObject(mod, "MassAction", (PyObject *) &MassActionType) < 0) {
    Py_DECREF(&MassActionType);
    Py_DECREF(mod);
    return NULL;
  }

  return mod;
}
