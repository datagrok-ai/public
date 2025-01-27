import { useEffect, useState, useTransition } from 'react';

async function readAllowedOrigins() {
  return [];
}

export function useWhiteList() {
  const [allowedOrigins, setAllowedOrigins] = useState<Array<string>>([]);
  const [isFetchAllowedOriginsPending, startTransition] = useTransition();

  useEffect(() => {
    startTransition(() => {
      void readAllowedOrigins().then((whitelist) =>
        setAllowedOrigins(whitelist),
      );
    });
  }, []);

  return {
    allowedOrigins,
    isFetchAllowedOriginsPending,
  };
}
