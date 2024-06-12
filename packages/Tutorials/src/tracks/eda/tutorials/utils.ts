export async function waitForElementClick(element: HTMLElement | null): Promise<void> {
    return new Promise<void>((resolve) => {
      const intervalId = setInterval(() => {
        if (!element) return;
        
        clearInterval(intervalId);
        element.addEventListener('click', () => resolve());
      }, 500);
    });
  }