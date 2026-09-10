// module-bar.js — module bar button -> nav pane switching
//
// Receives: { value: "gea" }
//
// Why a synthetic click instead of bslib::nav_select():
//   nav_select() -> shiny::updateTabsetPanel() -> shiny's BootstrapTabInputBinding
//   .setValue(), which on Bootstrap 5 (isBS3() === !window.bootstrap, and bslib
//   ships bootstrap.bundle.min.js so window.bootstrap is defined) looks for
//   `.nav-link` anchors. bslib emits legacy BS3 nav markup — `li > a` with NO
//   nav-link class, styled via `@extend .nav-link` in bs3compat/_nav_compat.scss
//   — so that selector matches nothing and setValue() silently no-ops.
//
//   bs3compat.js instead registers a delegated click handler on
//   '[data-bs-toggle="tab"]' that calls $(this).tab("show"), matching the markup
//   bslib actually produces. Clicking the anchor therefore goes through bslib's
//   own supported path. It works while the anchor is visually hidden: a
//   dispatched click bubbles to the document handler regardless of layout.
Shiny.addCustomMessageHandler("module_bar_select", function(msg) {
    var nav = document.getElementById("main_tabs");
    if (!nav || !msg || !msg.value) return;

    var anchor = nav.querySelector('a[data-value="' + msg.value + '"]');
    if (anchor) anchor.click();
});
