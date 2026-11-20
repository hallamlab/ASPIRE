document.addEventListener("DOMContentLoaded", () => {
    for (const block of document.querySelectorAll(".highlight-bash, .highlight-shell")) {
        const code = block.querySelector("pre");
        if (!code) continue;

        const button = document.createElement("button");
        button.type = "button";
        button.className = "copy-commands";
        button.title = "Copy commands";
        button.setAttribute("aria-label", "Copy commands");
        let resetTimer;
        button.addEventListener("click", async () => {
            clearTimeout(resetTimer);
            try {
                await navigator.clipboard.writeText(code.textContent.trimEnd());
                button.classList.add("copied");
                button.title = "Copied";
            } catch {
                const selection = window.getSelection();
                const range = document.createRange();
                range.selectNodeContents(code);
                selection.removeAllRanges();
                selection.addRange(range);
                button.title = "Commands selected. Copy manually.";
            }
            button.setAttribute("aria-label", button.title);
            resetTimer = setTimeout(() => {
                button.classList.remove("copied");
                button.title = "Copy commands";
                button.setAttribute("aria-label", button.title);
            }, 1800);
        });
        block.prepend(button);
    }
});
