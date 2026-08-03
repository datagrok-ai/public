$('#menubar').toggle();
$('#header-container').toggle();

define(['base/js/namespace', 'base/js/promises'], function(Jupyter, promises) {
  promises.notebook_loaded.then(function() {
    Jupyter.notebook.set_autosave_interval(3000);
  });
});

define(
  ['base/js/namespace', 'jquery'],
  function(jupyter, $) {
    $(jupyter.events).on("kernel_ready.Kernel", function () {
      let autorun = new URLSearchParams(window.location.search.slice(1)).get("autorun");
      if (autorun !== undefined && autorun !== null && autorun === 'true') {
        jupyter.actions.call('jupyter-notebook:run-all-cells-below');
        jupyter.actions.call('jupyter-notebook:save-notebook');
      }
    });
  }
);
