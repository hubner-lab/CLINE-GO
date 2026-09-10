// config-dirty.js — Shiny message handler for sidebar dirty-state badge
//
// Receives: { sidebar_id: "config_association", dirty_count: 3 }
// Toggles a .config-dirty-badge span inside the sidebar title element.
// The sidebar title selector targets the bslib sidebar panel's header.

// bio_chips_reset — sync chip active states when Reset is clicked
// Receives: { container_id, input_id, value: "bio_1,bio_3,bio_12" }
Shiny.addCustomMessageHandler("bio_chips_reset", function(msg) {
    var container = document.getElementById(msg.container_id);
    var inputEl   = document.getElementById(msg.input_id);
    if (!container) return;

    var selected = msg.value ? msg.value.split(",").map(function(s) { return s.trim(); }) : [];
    // If empty, default to all selected
    if (selected.length === 0) {
        for (var i = 1; i <= 19; i++) selected.push("bio_" + i);
    }

    container.querySelectorAll(".bc-chip").forEach(function(chip) {
        if (chip.dataset.invariant === "true") return;
        var active = selected.indexOf(chip.dataset.bio) !== -1;
        chip.classList.toggle("bc-active", active);
    });

    // Sync the hidden bridge input too
    if (inputEl) {
        inputEl.value = selected.join(",");
        Shiny.setInputValue(msg.input_id, selected.join(","), {priority: "event"});
    }
});

// pheno_chips_reset — sync phenotype trait chip active states when Reset is clicked
// Receives: { container_id, input_id, value: "height,flowering_time" }
Shiny.addCustomMessageHandler("pheno_chips_reset", function(msg) {
    var container = document.getElementById(msg.container_id);
    var inputEl   = document.getElementById(msg.input_id);
    if (!container) return;

    // Parse selected traits; empty = all selected
    var selected = msg.value ? msg.value.split(",").map(function(s) { return s.trim(); }).filter(Boolean) : [];
    var allChips = container.querySelectorAll(".bc-chip");
    if (selected.length === 0) {
        // Default: select all
        allChips.forEach(function(chip) { chip.classList.add("bc-active"); });
        var all = Array.from(allChips).map(function(c) { return c.dataset.trait; });
        if (inputEl) {
            inputEl.value = all.join(",");
            Shiny.setInputValue(msg.input_id, all.join(","), {priority: "event"});
        }
    } else {
        allChips.forEach(function(chip) {
            chip.classList.toggle("bc-active", selected.indexOf(chip.dataset.trait) !== -1);
        });
        if (inputEl) {
            inputEl.value = selected.join(",");
            Shiny.setInputValue(msg.input_id, selected.join(","), {priority: "event"});
        }
    }
});

// plot_rearm — re-arm the Manhattan loading overlay when the project changes.
// Receives: { id: "<ns>overlay" }
// Removes .plot-revealed so the overlay re-covers the new cold-start blink;
// the debounced plotly_afterplot listener baked into the widget re-reveals it.
Shiny.addCustomMessageHandler("plot_rearm", function(msg) {
    var out  = document.getElementById(msg.id);
    var wrap = out && out.closest(".plot-loading-wrap");
    if (!wrap) return;
    wrap.classList.remove("plot-revealed");
    delete wrap.dataset.revealFallback;   // re-arm max-timeout fallback
});

Shiny.addCustomMessageHandler("config_dirty_update", function(msg) {
    var sid  = msg.sidebar_id;
    var n    = msg.dirty_count;

    // The bslib sidebar toggle button has aria-controls = "{id}"
    // The title text span is inside the sidebar panel itself
    // Target: .bslib-sidebar-layout > .sidebar > .sidebar-title
    // Fall back to searching by sidebar id attribute
    var container = document.querySelector(
        '[data-bslib-sidebar-id="' + sid + '"], #' + sid
    );

    // Also find the sidebar toggle button (for badge on the button)
    var toggleBtn = document.querySelector('[aria-controls="' + sid + '"]');

    if (!container && !toggleBtn) return;

    // Manage badge on the toggle button
    if (toggleBtn) {
        var badge = toggleBtn.querySelector('.config-dirty-badge');
        if (n > 0) {
            if (!badge) {
                badge = document.createElement('span');
                badge.className = 'config-dirty-badge';
                toggleBtn.appendChild(badge);
            }
            badge.textContent = 'modified';
        } else {
            if (badge) badge.remove();
        }
    }

    // Also update the sidebar panel title if accessible
    if (container) {
        var titleEl = container.querySelector('.config-sidebar-title');
        if (titleEl) {
            var panelBadge = titleEl.querySelector('.config-dirty-badge');
            if (n > 0) {
                if (!panelBadge) {
                    panelBadge = document.createElement('span');
                    panelBadge.className = 'config-dirty-badge ms-2';
                    titleEl.appendChild(panelBadge);
                }
                panelBadge.textContent = 'modified';
            } else {
                if (panelBadge) panelBadge.remove();
            }
        }
    }
});
